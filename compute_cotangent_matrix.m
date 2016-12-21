function L = compute_cotangent_matrix(vertices,faces,varargin)
% COTANGENT compute the cotangents of each angle in mesh (V,F), more details
  % can be found in Section 1.1 of "Algorithms and Interfaces for Real-Time
  % Deformation of 2D and 3D shapes" [Jacobson 2013]
  % 
  % C = cotangent(V,F)
  % C = cotangent(V,F,'ParameterName',parameter_value,...)
  %
  % Known bugs:
  %   This seems to return 0.5*C and for tets already multiplies by
  %   edge-lengths
  %
  % Inputs:
  %   V  #V by dim list of rest domain positions
  %   F  #F by {3|4} list of {triangle|tetrahedra} indices into V
  %   Optional (3-manifolds only):
  %     'SideLengths' followed by #F by 3 list of edge lengths corresponding
  %       to: 23 31 12. In this case V is ignored.
  %       or
  %       followed by #T by 6 list of tet edges lengths:
  %       41 42 43 23 31 12
  %     'FaceAreas' followed by #T by 4 list of tet face areas
  % Outputs:
  %   C  #F by {3|6} list of cotangents corresponding
  %     angles for triangles, columns correspond to edges 23,31,12
  %     dihedral angles *times opposite edge length* over 6 for tets, 
  %       WRONG: columns correspond to edges 23,31,12,41,42,43 
  %       RIGHT: columns correspond to *faces* 23,31,12,41,42,43
  %
  % See also: cotmatrix
  %
  % Copyright 2013, Alec Jacobson (jacobson@inf.ethz.ch)
  %


% default values
l = [];
% triangles
if isempty(l)
  % edge lengths numbered same as opposite vertices
  l = [ ...
    sqrt(sum((vertices(faces(:,2),:)-vertices(faces(:,3),:)).^2,2)) ...
    sqrt(sum((vertices(faces(:,3),:)-vertices(faces(:,1),:)).^2,2)) ...
    sqrt(sum((vertices(faces(:,1),:)-vertices(faces(:,2),:)).^2,2)) ...
    ];
end

l1 = l(:,1); l2 = l(:,2); l3 = l(:,3);
% semiperimeters
s = (l1 + l2 + l3)*0.5;
% Heron's formula for area
dblA = 2*sqrt( s.*(s-l1).*(s-l2).*(s-l3));
% cotangents and diagonal entries for element matrices
% correctly divided by 4 (alec 2010)
Cot = [(l2.^2 + l3.^2 -l1.^2)./dblA/4 ...
     (l1.^2 + l3.^2 -l2.^2)./dblA/4 ...
     (l1.^2 + l2.^2 -l3.^2)./dblA/4 ...
    ];

L = sparse(faces(:,[2 3 1]), faces(:,[3 1 2]), Cot,size(vertices,1),size(vertices,1));
L = L+L';
L = L-diag(sum(L,2));
