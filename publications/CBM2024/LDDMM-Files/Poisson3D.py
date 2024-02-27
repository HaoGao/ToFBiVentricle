import numpy as np
import torch
import math
from torch.autograd import Variable
from numpy import loadtxt
from scipy.sparse import csr_matrix
from scipy.sparse import spdiags
from scipy.sparse.linalg import lsqr


from core import default



class Poisson3D:
    ####################################################################################################################
    ### Constructor:
    ####################################################################################################################

    def __init__(self, node,elem,boundary_node):
        self.points = node
        self.source = elem
        self.bzd = boundary_node

        ####################################################################################################################
        ### This part is to include Possion Function to construct 3D Four-node tetrahedral element
        ### in solid geometry
        #################################################################################################################### 
    def Poisson3(node,elem,boundary_node):


        def shapefunction(node,elem):

            NT = elem.shape[0] 
            face=np.concatenate((elem[:,[1, 3, 2]],  elem[:,[0, 2, 3]], \
                 elem[:,[0, 3, 1]], elem[:,[0, 1, 2]]), axis=0)
            v12 = node[face[:,1],:]-node[face[:,0],:]
            v13 = node[face[:,2],:]-node[face[:,0],:]
            normal = np.cross(v12,v13)
            v12 = v12[3*NT:4*NT,:]
            v13 = v13[3*NT:4*NT,:]
            v14 = node[elem[:,3],:]-node[elem[:,0],:]
            volume = np.sum(np.cross(v12,v13)*v14, axis=1)/6.0  
            volume=np.absolute(volume)
            Dlambda = np.zeros((NT,3,4))
            volume0=np.column_stack((6.0*volume, 6.0*volume,6.0*volume))
            Dlambda[0:NT,:,0] = normal[0:NT,:]/volume0    
            Dlambda[0:NT,:,1] = normal[NT:2*NT,:]/volume0    
            Dlambda[0:NT,:,2] = normal[2*NT:3*NT,:]/volume0    
            Dlambda[0:NT,:,3] = normal[3*NT:4*NT,:]/volume0    

            return Dlambda, volume 


        #@staticmethod
        def getboundary(b,boundary_node,A,Ndof):

            u = np.zeros((Ndof,1))
            fixedNode = [] 
            freeNode = []

            ## Part 1: Modify the matrix for Dirichlet condition
            # Find Dirichlet boundary nodes: fixedNode
            for i in range(boundary_node.shape[0]):
                if boundary_node[i,0]==1:
                    fixedNode.append(i)
                else:
                    freeNode.append(i)

            # Modify the matrix for different boundary conditions   
            # Dirichlet boundary condition
    
        
            bdidx = np.zeros((Ndof,1))
            bdidx[fixedNode] = 1
            Tbd = spdiags(bdidx.T,0,Ndof,Ndof)
            T = spdiags((1-bdidx).T,0,Ndof,Ndof)
            AD = T*A*T + Tbd

            # Dirichlet boundary condition
            g_D=boundary_node[:,1]
            g_D=g_D.reshape(len(g_D),1)
            u[fixedNode] =g_D[fixedNode]
    
            b = b - A*u
            b[fixedNode] = u[fixedNode]

            return AD,b,u,freeNode


        # Set up optional input arguments  
        N = node.shape[0] 
        NT = elem.shape[0]
        Ndof = N
        # Compute geometric quantities and gradient of local basis
        Dphi,volume = shapefunction(node,elem)
        # Assemble stiffness matrix
        A = csr_matrix((Ndof,Ndof), dtype=np.float16)
        for i in range(4):
            for j in range(i,4):
                Aij = np.sum(Dphi[:,:,i]*Dphi[:,:,j], axis=1)*volume
                if j==i:
                    A = A + csr_matrix((Aij, (elem[:,i], elem[:,j])), shape=(Ndof, Ndof))
                else:
                    AijN=np.concatenate((Aij,Aij), axis=0)
                    elemN1=np.concatenate((elem[:,i],elem[:,j]), axis=0)
                    elemN2=np.concatenate((elem[:,j],elem[:,i]), axis=0)
                    A = A + csr_matrix((AijN, (elemN1, elemN2)), shape=(Ndof, Ndof))
        # Assemble right hand side
        b = np.zeros((Ndof,1))

        # Set up boundary conditions
        AD,b,u,freeNode = getboundary(b,boundary_node,A,Ndof)

        # Solve the system of linear equations
        # solve
        xe=AD[freeNode]
        xd=xe[:,freeNode]
        x=lsqr(xd,b[freeNode])[0]
        x=x.reshape(len(x),1)
        u[freeNode]=x

       # Compute Du
        dudx = np.multiply(u[elem[:,0]],Dphi[:,0,0].reshape(len(elem[:,0]),1))\
               +np.multiply(u[elem[:,1]],Dphi[:,0,1].reshape(len(elem[:,0]),1))\
               +np.multiply(u[elem[:,2]],Dphi[:,0,2].reshape(len(elem[:,0]),1)) \
               +np.multiply(u[elem[:,3]],Dphi[:,0,3].reshape(len(elem[:,0]),1)) 
        dudy = np.multiply(u[elem[:,0]],Dphi[:,1,0].reshape(len(elem[:,0]),1))\
               +np.multiply(u[elem[:,1]],Dphi[:,1,1].reshape(len(elem[:,0]),1)) \
               +np.multiply(u[elem[:,2]],Dphi[:,1,2].reshape(len(elem[:,0]),1))\
               +np.multiply(u[elem[:,3]],Dphi[:,1,3].reshape(len(elem[:,0]),1)) 
        dudz = np.multiply(u[elem[:,0]],Dphi[:,2,0].reshape(len(elem[:,0]),1))\
               +np.multiply(u[elem[:,1]],Dphi[:,2,1].reshape(len(elem[:,0]),1)) \
               +np.multiply(u[elem[:,2]],Dphi[:,2,2].reshape(len(elem[:,0]),1))\
               +np.multiply(u[elem[:,3]],Dphi[:,2,3].reshape(len(elem[:,0]),1)) 
        Du = np.concatenate((dudx, dudy, dudz), axis=1)

        return u, Du
    
  



