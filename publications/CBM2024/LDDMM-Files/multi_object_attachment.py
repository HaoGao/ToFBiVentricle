import numpy as np
import torch
import math
from torch.autograd import Variable
from numpy import loadtxt
from scipy.sparse import csr_matrix
from scipy.sparse import spdiags
from scipy.sparse.linalg import lsqr

from core import default
from core.model_tools.attachments.ModificationFromDBG import ModificationFromDBG


class MultiObjectAttachment:
    ####################################################################################################################
    ### Constructor:
    ####################################################################################################################

    def __init__(self, attachment_types, kernels):
        # List of strings, e.g. 'varifold' or 'current'.
        self.attachment_types = attachment_types
        # List of kernel objects.
        self.kernels = kernels


    ####################################################################################################################
    ### Public methods:
    ####################################################################################################################

    def compute_weighted_distance(self, data, multi_obj1, multi_obj2, inverse_weights):
        """
        Takes two multiobjects and their new point positions to compute the distances
        """
        distances = self.compute_distances(data, multi_obj1, multi_obj2)
        assert distances.size()[0] == len(inverse_weights)
        inverse_weights_torch = torch.from_numpy(np.array(inverse_weights)).type(list(data.values())[0].type())
        return torch.sum(distances / inverse_weights_torch)

    def compute_distances(self, data, multi_obj1, multi_obj2):
        """
        Takes two multiobjects and their new point positions to compute the distances.
        """
        assert len(multi_obj1.object_list) == len(multi_obj2.object_list), \
            "Cannot compute distance between multi-objects which have different number of objects"
        distances = torch.zeros((len(multi_obj1.object_list),)).type(list(data.values())[0].type())

        pos = 0
        for i, obj1 in enumerate(multi_obj1.object_list):
            obj2 = multi_obj2.object_list[i]

            if self.attachment_types[i].lower() == 'current':
                distances[i] = self.current_distance(
                    data['landmark_points'][pos:pos + obj1.get_number_of_points()], obj1, obj2, self.kernels[i])
                pos += obj1.get_number_of_points()

            elif self.attachment_types[i].lower() == 'pointcloud':
                distances[i] = self.point_cloud_distance(
                    data['landmark_points'][pos:pos + obj1.get_number_of_points()], obj1, obj2, self.kernels[i])
                pos += obj1.get_number_of_points()

            elif self.attachment_types[i].lower() == 'varifold':
                distances[i] = self.varifold_distance(
                    data['landmark_points'][pos:pos + obj1.get_number_of_points()], obj1, obj2, self.kernels[i])
                pos += obj1.get_number_of_points()

            elif self.attachment_types[i].lower() == 'landmark':
                distances[i] = self.landmark_distance(
                    data['landmark_points'][pos:pos + obj1.get_number_of_points()], obj2)
                pos += obj1.get_number_of_points()

            elif self.attachment_types[i].lower() == 'l2':
                assert obj1.type.lower() == 'image' and obj2.type.lower() == 'image'
                distances[i] = self.L2_distance(data['image_intensities'], obj2)

            else:
                assert False, "Please implement the distance {e} you are trying to use :)".format(
                    e=self.attachment_types[i])

        return distances

    ####################################################################################################################
    ### Auxiliary methods:
    ####################################################################################################################

    @staticmethod
    def current_distance(points, source, target, kernel):
        """
        Compute the current distance between source and target, assuming points are the new points of the source
        We assume here that the target never moves.
        """

        c1, n1, c2, n2 = MultiObjectAttachment.__get_source_and_target_centers_and_normals(points, source, target)

        def current_scalar_product(points_1, points_2, normals_1, normals_2):
            return torch.dot(normals_1.view(-1), kernel.convolve(points_1, points_2, normals_2).view(-1))

        if target.norm is None:
            target.norm = current_scalar_product(c2, c2, n2, n2)

        return current_scalar_product(c1, c1, n1, n1) + target.norm - 2 * current_scalar_product(c1, c2, n1, n2)

    @staticmethod
    def point_cloud_distance(points, source, target, kernel):
        """
        Compute the point cloud distance between source and target, assuming points are the new points of the source
        We assume here that the target never moves.
        """

        c1, n1, c2, n2 = MultiObjectAttachment.__get_source_and_target_centers_and_normals(points, source, target)

        def point_cloud_scalar_product(points_1, points_2, normals_1, normals_2):
            return torch.dot(normals_1.view(-1),
                             kernel.convolve(points_1, points_2, normals_2, mode='pointcloud').view(-1))

        if target.norm is None:
            target.norm = point_cloud_scalar_product(c2, c2, n2, n2)

        return point_cloud_scalar_product(c1, c1, n1, n1) + target.norm - 2 * point_cloud_scalar_product(c1, c2, n1, n2)

    @staticmethod
    def varifold_distance(points, source, target, kernel):

        """
        Returns the current distance between the 3D meshes
        source and target are SurfaceMesh objects
        points are source points (torch tensor)
        """

        c1, n1, c2, n2 = MultiObjectAttachment.__get_source_and_target_centers_and_normals(points, source, target)
        

        # alpha = normales non unitaires
        areaa = torch.norm(n1, 2, 1)
        areab = torch.norm(n2, 2, 1)

        nalpha = n1 / areaa.unsqueeze(1)
        nbeta = n2 / areab.unsqueeze(1)

        #print (source.connectivity)

        def varifold_scalar_product(x, y, areaa, areab, nalpha, nbeta):
            return torch.dot(areaa.view(-1), kernel.convolve((x, nalpha), (y, nbeta), areab.view(-1, 1),
                                                             mode='varifold').view(-1))

        if target.norm is None:
            target.norm = varifold_scalar_product(c2, c2, areab, areab, nbeta, nbeta)
        
        ####################################################################
        #***** Modification DGB 2022/06/08
        bzd=varifold_scalar_product(c1, c1, areaa, areaa, nalpha, nalpha) + target.norm \
            - 2 * varifold_scalar_product(c1, c2, areaa, areab, nalpha, nbeta) 
        
        ####################################################################
        #***** Modification DGB 2022/06/08, edge length 
        #sum_variance_length=ModificationFromDBG.edge_lenght_variance_dbg(points, source, bzd)
        #sum_variance_length=0.0
        ####################################################################
        #***** Modification DGB 2022/06/12, Possion function, element quality
        #sum_variance_quality_element=ModificationFromDBG.quality_solid_element(points, source, bzd)
        ####################################################################
        #***** Modification DGB 2022/06/12, Possion function, element quality
        F_undia, quality_variance, volume_variance=ModificationFromDBG.quality_axial_growth(points, source, bzd)
        #quality_variance=0.0

        print('F_undia',F_undia)
        print('quality_variance', quality_variance)   
        print('volume_variance', volume_variance)      
        print('Vold_distance', bzd)
        print('total', bzd+F_undia+quality_variance+volume_variance)
        #***** Modification DGB 2022/06/08
        return bzd+F_undia+quality_variance+volume_variance

    @staticmethod
    def landmark_distance(points, target):
        """
        Point correspondance distance
        """
        target_points = target.get_points_torch(tensor_scalar_type=points.type())
        return torch.sum((points.view(-1) - target_points.view(-1)) ** 2)

    @staticmethod
    def L2_distance(intensities, target):
        """
        L2 image distance.
        """
        target_intensities = target.get_intensities_torch(tensor_scalar_type=intensities.type())
        return torch.sum((intensities.view(-1) - target_intensities.view(-1)) ** 2)

    ####################################################################################################################
    ### Private methods:
    ####################################################################################################################

    @staticmethod
    def __get_source_and_target_centers_and_normals(points, source, target):
        tensor_scalar_type = points.type()
        tensor_integer_type = {
            'cpu': 'torch.LongTensor',
            'cuda': 'torch.cuda.LongTensor'
        }[points.device.type]

        c1, n1 = source.get_centers_and_normals(points,
                                                tensor_scalar_type=tensor_scalar_type,
                                                tensor_integer_type=tensor_integer_type)
        c2, n2 = target.get_centers_and_normals(tensor_scalar_type=tensor_scalar_type,
                                                tensor_integer_type=tensor_integer_type)
      
        return c1, n1, c2, n2
        

    
       




