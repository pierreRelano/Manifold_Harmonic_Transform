path(path, 'mesh/');
path(path, 'mesh/meshes/');
path(path, 'mesh/toscahires-mat/');
path(path, 'mesh/mesh_deform/');
path(path, 'mesh/Models/');
path(path, 'lib/gptoolbox/mesh');
path(path, 'lib/gptoolbox/external/toolbox_fast_marching/toolbox');
path(path, 'lib/toolbox_graph/');% Mesh Smoothing  with Manifold Harmonics [Vallet et al. 2008]


clc
clear all
close all


% READ MESH FILE
filePath = 'mesh/Geometry/icosahedron.off';
[~, remain] = strtok(filePath, '.');
if strcmp('.obj', remain) == 1
     [V, F] = readOBJ(filePath);     %Read OBJ file
elseif strcmp('.off', remain) == 1
    [V, F] = readOFF(filePath);  % Read OFF file
else
    disp('Unknown file type!');
end

% renaming indices of vertices of triangles for convenience
i1 = F(:,1); 
i2 = F(:,2); 
i3 = F(:,3); 
% #F x 3 matrices of triangle edge vectors, named after opposite vertices
v1 = V(i3,:) - V(i2,:);  %edges vector coordinates opposed to the angles of interest
v2 = V(i1,:) - V(i3,:); 
v3 = V(i2,:) - V(i1,:);
% computing *unsigned* areas 
if size(V,2) == 2
    % 2d vertex data
    dblA = abs(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1));
elseif size(V,2) == 3
    %n  = cross(v1,v2,2);  dblA  = multinorm(n,2);

    % area of parallelogram is twice area of triangle
    % area of parallelogram is || v1 x v2 || 
    n  = cross(v1,v2,2); 
    % THIS DOES MATRIX NORM!!! don't use it!!
    % dblA  = norm(n,2);

    % This does correct l2 norm of rows
    dblA = (sqrt(sum((n').^2)))';
else 
    error('unsupported vertex dimension %d', size(V,2))
end
% cotangents and diagonal entries for element matrices
cot12 = -dot(v1,v2,2)./dblA/2; 
cot23 = -dot(v2,v3,2)./dblA/2; 
cot31 = -dot(v3,v1,2)./dblA/2;
% diag entries computed from the condition that rows of the matrix sum up to 1
% (follows from  the element matrix formula E_{ij} = (v_i dot v_j)/4/A )
diag1 = -cot12-cot31; 
diag2 = -cot12-cot23; 
diag3 = -cot31-cot23;
% indices of nonzero elements in the matrix for sparse() constructor
i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
% values corresponding to pairs form (i,j)
v = [cot12 cot12 cot23 cot23 cot31 cot31 diag1 diag2 diag3];
% for repeated indices (i,j) sparse automatically sums up elements, as we
% want
L = full(sparse(i,j,v,size(V,1),size(V,1)));