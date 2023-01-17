%% a possion 2D example using FE solver
clear all; close all; clc;
path(path, '../../src');

%% load the mesh 
nodes_matrix = load('nodes.txt');
nodes_matrix = nodes_matrix(:, 1:2);
elems_matrix = load('elements.txt');
elems_matrix = elems_matrix(:, 4:6) + 1; % index starts at 1


Ndof = size(nodes_matrix,1); % total degree of freedom
fxy = 1;
fbc = 1;
epsilon = 1e-6;
penalty = 1e6;

% stiffness matrix and the right hand side
K = zeros([Ndof, Ndof]);
f = zeros([Ndof, 1]);

% loop all elements
[gpcoor, gpweight]= gp_tri_lin;
gpNumber = length(gpweight);
for el = 1 : size(elems_matrix, 1)
    
    elem = elems_matrix(el,:);
    elem_dofs = dof_map_element(elem);
    n_dofs = length(elem_dofs);
    
    XCOOR = nodes_matrix(elem,:);
    
    Ktemp = zeros([n_dofs, n_dofs]);
    ftemp = zeros([n_dofs,1]);
    for qp = 1 : gpNumber
        [shape, dshape, Jxw] = shape_tri_lin(XCOOR,gpcoor(qp,:));
        
        %% assemble the K matrix
        for i = 1 : n_dofs
            dphi_i = dshape(i,:);
            for j = 1 : n_dofs
                dphi_j = dshape(j,:);
                Ktemp(i,j) = Ktemp(i,j) + Jxw*dot(dphi_i, dphi_j);
            end % j
        end %i
        
        %% assemble the right hand side
        for i = 1 : n_dofs
            phi_i = shape(i);
            ftemp(i) = ftemp(i) +  Jxw*phi_i*(-fxy);
        end
        
    end % end for qp
    
    K = constrain_stiffness_matrix(K, Ktemp, elem_dofs);
    f = constrain_rhs(f, ftemp, elem_dofs);
    
end% end for element


%% set up the boundary conditions based on nodes
for nl = 1 : size(nodes_matrix, 1)
    X = nodes_matrix(nl,:);
    if abs(X(1)-1.0)<epsilon || abs(X(1)+1.0)<epsilon ... 
       || abs(X(2)-1.0)<epsilon || abs(X(2)+1.0)<epsilon
        %at the boundary
        ndof = dof_map_node(nl);
        K(ndof, ndof) = K(ndof, ndof) + penalty;
        f(ndof)  = f(ndof) + penalty*fbc;
    end
end

% solve the system
u = inv(K)*f;

%% plot the result
h3d = figure(); hold on;
patch('Faces',elems_matrix,'Vertices',nodes_matrix, 'FaceVertexCData',u,...
      'FaceColor','interp', 'EdgeColor', 'none');
colormap(h3d,jet);
colorbar;
axis equal;

