import numpy as np
import torch
import math
from torch.autograd import Variable
from numpy import loadtxt
from scipy.sparse import csr_matrix
from scipy.sparse import spdiags
from scipy.sparse.linalg import lsqr

from core import default
from core.model_tools.attachments.Poisson3D import Poisson3D




class ModificationFromDBG:
    ####################################################################################################################
    ### Constructor:
    ####################################################################################################################

    def __init__(self, points, source, bzd):
        self.points = points
        self.source = source
        self.bzd = bzd

    
        ###################################################################################

        #***** Modification DGB 2022/09/16
        ####################################################################################################################
        ### Add function to compute growth amount along three axial direction
        ####################################################################################################################
   
    def quality_axial_growth(points, source, bzd): 
        #read mesh data and boundary conditions
        elem=source.connectivity 
        node_intial=source.points  #surface nodes
        node_updated=points.detach().numpy()        #surface nodes
    
        bd_node_surface=node_updated-node_intial #surface nodes displacement
        # Map surface nodes and elemnt back to the solid nodes 
        # should only read once at the begining, later do this 
        file_path='/home/pgrad2/2306902g/DBGuan/DB/DTMRImapping/RatHeart_Growth_MCT/data/'
        #file_path='C:/Users/dg215r/OneDrive - University of Glasgow/NoWay/LianTian/'+\
         #         'Heart modeling/HeartRat/Python_poisson/AxialGrowth/'
        #bd_node_surface=loadtxt(file_path+"Template_Target_dxdydz.txt", delimiter="\t", unpack=False) 
        bd_node_map=loadtxt(file_path+"Template_surface_solid_node_map.txt", delimiter="\t", unpack=False) 
        solid_elem=loadtxt(file_path+"TemplateElement.txt", delimiter="\t", unpack=False) 
        solid_node=loadtxt(file_path+"TemplateNode.txt", delimiter="\t", unpack=False) 
        solid_elem_LV_RV_SEP=loadtxt(file_path+"Template_LVRVSEP.txt", delimiter="\t", unpack=False) 
        fsn_sys=loadtxt(file_path+"274_fibre.txt", delimiter=",", unpack=False) 
        sheetboundary=loadtxt(file_path+"SheetBoundary.txt", delimiter="\t", unpack=False)
        solid_elem=solid_elem.astype(int)-1
        bd_node_map=bd_node_map.astype(int)
        solid_elem_LV_RV_SEP=solid_elem_LV_RV_SEP.astype(int)
  

        Nnode=bd_node_map.shape[0]
        Nelem=solid_elem.shape[0]
        
        bd_node_solid_x=np.zeros((Nnode,2))
        bd_node_solid_y=np.zeros((Nnode,2))
        bd_node_solid_z=np.zeros((Nnode,2))

        for ia, ib, ic in bd_node_map:
            if ia==1:
                bd_node_solid_x[ib-1,0]=1
                bd_node_solid_x[ib-1,1]=bd_node_surface[ic-1,0] 
                bd_node_solid_y[ib-1,0]=1
                bd_node_solid_y[ib-1,1]=bd_node_surface[ic-1,1]   
                bd_node_solid_z[ib-1,0]=1
                bd_node_solid_z[ib-1,1]=bd_node_surface[ic-1,2]  
                
        # call poisson solver
                #x
        soln,Du = Poisson3D.Poisson3(solid_node,solid_elem,bd_node_solid_x)
        u_x=soln.reshape(len(soln),1)
                #y
        soln,Du = Poisson3D.Poisson3(solid_node,solid_elem,bd_node_solid_y)
        u_y=soln.reshape(len(soln),1)
                #z
        soln,Du = Poisson3D.Poisson3(solid_node,solid_elem,bd_node_solid_z)
        u_z=soln.reshape(len(soln),1)
                #
        u_xyz=np.concatenate((u_x,u_y,u_z), axis=1)
        
        ######### reconstruct solid geometry with updated node positions

        node_deform=solid_node+u_xyz

        soln,Du = Poisson3D.Poisson3(node_deform,solid_elem,sheetboundary)
        
        
        ###################
        ######### compute deformation gradient tensor 
        def IsoTet4ShapeFunDer(xyztet):

            x1=xyztet[0,0]
            x2=xyztet[0,1]
            x3=xyztet[0,2]
            x4=xyztet[0,3]
            y1=xyztet[1,0]
            y2=xyztet[1,1]
            y3=xyztet[1,2]
            y4=xyztet[1,3]
            z1=xyztet[2,0]
            z2=xyztet[2,1]
            z3=xyztet[2,2]
            z4=xyztet[2,3]

            x12=x1-x2
            x13=x1-x3
            x14=x1-x4
            x23=x2-x3
            x24=x2-x4
            x34=x3-x4
            x21=-x12
            x31=-x13
            x41=-x14
            x32=-x23
            x42=-x24
            x43=-x34
            y12=y1-y2
            y13=y1-y3
            y14=y1-y4
            y23=y2-y3
            y24=y2-y4
            y34=y3-y4
            y21=-y12
            y31=-y13
            y41=-y14
            y32=-y23
            y42=-y24
            y43=-y34
            z12=z1-z2
            z13=z1-z3
            z14=z1-z4
            z23=z2-z3
            z24=z2-z4
            z34=z3-z4
            z21=-z12
            z31=-z13
            z41=-z14
            z32=-z23
            z42=-z24
            z43=-z34

            #compute the a b c
            abc=np.zeros((4,3))
            abc[0,0]=y42*z32-y32*z42
            abc[0,1]=x32*z42-x42*z32
            abc[0,2]=x42*y32-x32*y42
            abc[1,0]=y31*z43-y34*z13
            abc[1,1]=x43*z31-x13*z34
            abc[1,2]=x31*y43-x34*y13
            abc[2,0]=y24*z14-y14*z24
            abc[2,1]=x14*z24-x24*z14
            abc[2,2]=x24*y14-x14*y24
            abc[3,0]=y13*z21-y12*z31
            abc[3,1]=x21*z13-x31*z12
            abc[3,2]=x13*y21-x12*y31

            # compute the volume
            V01=x2*(y3*z4-y4*z3)+x3*(y4*z2-y2*z4)+x4*(y2*z3-y3*z2);
            V02=x1*(y4*z3-y3*z4)+x3*(y1*z4-y4*z1)+x4*(y3*z1-y1*z3);
            V03=x1*(y2*z4-y4*z2)+x2*(y4*z1-y1*z4)+x4*(y1*z2-y2*z1);
            V04=x1*(y3*z2-y2*z3)+x2*(y1*z3-y3*z1)+x3*(y2*z1-y1*z2);
            Vcol=(V01+V02+V03+V04)/6.0;    
            
            return abc, Vcol
          
        


        ###################################################################
        ## element quality_
        ###################################################################

        solid_elem_new=np.concatenate((solid_elem,solid_elem[:,0:3]), axis=1)
        higness=np.zeros((Nelem,4))
        volumer=np.zeros((Nelem,4))
        
        for i in range(4):
            ag = torch.from_numpy(node_deform[solid_elem_new[:, i]])
            bg = torch.from_numpy(node_deform[solid_elem_new[:, i+1]])
            cg = torch.from_numpy(node_deform[solid_elem_new[:, i+2]])
            dg = torch.from_numpy(node_deform[solid_elem_new[:, i+3]])  
            length1=torch.norm(ag-bg,dim=1)
            length2=torch.norm(bg-cg,dim=1)
            length3=torch.norm(cg-ag,dim=1)
            normals = torch.cross(ag-bg, bg-cg)
            normals=torch.div(torch.t(normals), torch.norm(normals,dim=1))
            normals=torch.t(normals)
            le_mean=(length1+length2+length3)/2.
            area=torch.sqrt(le_mean*(le_mean-length1)*(le_mean-length2)*(le_mean-length3))
            #back to numpy
            gd=(dg-bg).detach().numpy()
            gt=normals.detach().numpy()
            agrea=area.detach().numpy()
            for tm in range(Nelem):
                higness[tm,i]=np.dot(gd[tm],gt[tm])/math.sqrt(agrea[tm])/1.24
                volumer[tm,i]=np.dot(gd[tm],gt[tm])*agrea[tm]/3.0
    
        higness=abs(higness)        
        higness_quality=abs(1.0-higness.min(axis=1)) #the samll the betteer


        ##################################################
        dxdydz=np.zeros((3,4))
        xyztet=np.zeros((3,4))
        #DeformationGrad=np.zeros((Nelem,3,3))
#######    non-diagonal components
        F_undia=0.0
        sigma=0.0
        F_se=0.0
#######    volume_variance
        
        volume_LVFW=[]
        volume_RVFW=[]
        volume_SEP=[]
#######    quality 
        higness_quality_sum=0.0
        diver=0.0
        v_note=0.0

        for i in range(Nelem):
            #for element quality
            if higness_quality[i] > 0.5:
                higness_quality_sum=higness_quality_sum+higness_quality[i]

################################
            xyztet=solid_node[solid_elem[i,:],:].T
            dxdydz=u_xyz[solid_elem[i,:],:].T

    
            abc, Vcol=IsoTet4ShapeFunDer(xyztet)
    
	        ###the deformation gradient tensor
            F=np.zeros((3,3))
            F=np.identity(3)+np.matmul(dxdydz, abc)/6.0/Vcol
###########################################
            ##compare gradient and sheet fibre
            sheet=np.matmul(F,fsn_sys[i,3:6])
            sheet=sheet/np.linalg.norm(sheet)
            gradi=Du[i,:]/np.linalg.norm(Du[i,:])
            diver=diver+1.0-abs(np.vdot(sheet,gradi))
            #print(sheet,gradi,diver)


###############    volume 
            volumes=np.linalg.det(F)     
            #print (volumes)

        #for i in range(Nelem):
            # for element volume
            if solid_elem_LV_RV_SEP[i]==1:  #LVFW
                volume_LVFW=np.append(volume_LVFW, volumes)
            elif solid_elem_LV_RV_SEP[i]==2:  #RVFW
                volume_RVFW=np.append(volume_RVFW, volumes)	
            else:
                volume_SEP=np.append(volume_SEP, volumes)

                
        volume_LVFW=abs(volume_LVFW)
        volume_RVFW=abs(volume_RVFW)
        volume_SEP=abs(volume_SEP)
        #print (len(volume_LVFW)+len(volume_RVFW)+len(volume_SEP))
        
        volume_LVFW_mean=np.mean(volume_LVFW)
        volume_RVFW_mean=np.mean(volume_RVFW)  
        volume_SEP_mean=np.mean(volume_SEP) 
        # Normalised to mean volume of each region
        variance_LVFW=np.square(volume_LVFW/volume_LVFW_mean-1.0)        
        variance_RVFW=np.square(volume_RVFW/volume_RVFW_mean-1.0)   
        variance_SEP=np.square(volume_SEP/volume_SEP_mean-1.0)  
                
        volume_variance=variance_LVFW.sum() +variance_RVFW.sum() +variance_SEP.sum()+v_note     
            
        #return F_undia

######################################################################
######################################################################

        scale_ratio=0.2 #1000.0*
        scale_ratio2=0.2
        #scale_ratio3=100.0 #03/04/2023 test for iteration
        quality_variance=scale_ratio2*10000.0*scale_ratio/Nelem*higness_quality_sum  #0.5 for each case, #0.2 for 17 regions #0.2 for 3 regions
        volume_variance=0.0*5000.0*scale_ratio/Nelem*volume_variance
        #F_undia=0.0*scale_ratio/Nelem*F_undia+sigma*scale_ratio
        #quality_variance=diver**scale_ratio*0.0
        F_undia=0.0 #sigma*scale_ratio
        #volume_variance=0.0
            
        return F_undia, quality_variance, volume_variance
        ###################################################################################

