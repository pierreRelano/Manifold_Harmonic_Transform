function L = compute_geometric_laplacian(vertex,face,options)

% compute_mesh_weight - compute a weight matrix
%
%   W = compute_mesh_weight(vertex,face,type,options);
%
%   W is sparse weight matrix and W(i,j)=0 is vertex i and vertex j are not
%   connected in the mesh.
%
%   type is either 
%       'combinatorial': W(i,j)=1 is vertex i is conntected to vertex j.
%       'distance': W(i,j) = 1/d_ij^2 where d_ij is distance between vertex
%           i and j.
%       'conformal': W(i,j) = cot(alpha_ij)+cot(beta_ij) where alpha_ij and
%           beta_ij are the adjacent angle to edge (i,j)
%
%   If options.normalize=1, the the rows of W are normalize to sum to 1.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if isfield(options, 'normalize')
    normalize = options.normalize;
else
    normalize = 1;
end
if isfield(options, 'symmetrize')
    symmetrize = options.symmetrize;
else
    symmetrize = 1;
end

options.normalize = 0;
[vertex,face] = check_face_vertex(vertex,face);
n = max(max(face));
verb = getoptions(options, 'verb', n>5000);


% geometric laplacian        
W = zeros(n,n);
ring = compute_vertex_face_ring(face);
for i = 1:n
    if verb
        progressbar(i,n);
    end
    for b = ring{i}
        % b is a face adjacent to a
        bf = face(:,b);
        % compute complementary vertices
        if bf(1)==i
            v = bf(2:3);
        elseif bf(2)==i
            v = bf([1 3]);
        elseif bf(3)==i
            v = bf(1:2);
        else
            error('Problem in face ring.');
        end
        j = v(1); k = v(2);
        vi = vertex(:,i);
        vj = vertex(:,j);
        vk = vertex(:,k);
        % angles
        alpha = compute_angle(vk-vi,vk-vj);
        beta = compute_angle(vj-vi,vj-vk);
        % add weight
        W(i,j) = W(i,j) + cot( alpha );
        W(i,k) = W(i,k) + cot( beta );
    end
end
if isfield(options, 'normalize') && options.normalize==1
    W = diag(sum(W,2).^(-1)) * W;
end

L = diag(sum(W,2)) - W;
% 
% if symmetrize==1 && normalize==0
%     L = diag(sum(W,2)) - W;
% elseif symmetrize==1 && normalize==1
%     L = speye(n) - diag(sum(W,2).^(-1/2)) * W * diag(sum(W,2).^(-1/2));
% elseif symmetrize==0 && normalize==1
%     L = speye(n) - diag(sum(W,2).^(-1)) * W;
% else
%     error('Does not work with symmetrize=0 and normalize=0');    
% end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha = compute_angle(u,v)
du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,eps); dv = max(dv,eps);
alpha = acos( sum(u.*v) / (du*dv) );
