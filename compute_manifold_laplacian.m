function L = compute_manifold_laplacian(vertex, faces, type, options)
%The combinatorial laplacian is a linear operator (thus a NxN matrix where N is the number of vertices). 
%It depends only on the connectivity of the mesh, thus on face only.

    if nargin<3
        type = 'contan';
    end

    %TODO: retest the indice of combinatorial laplacian
    switch lower(type)
        case 'combinatorial'
            %Compute edge list.
            n = size(vertex,2);
            E = [faces(:, [1 2]) faces(:, [2 3]) faces(:, [3 1])];
            p = size(E,2);
            %Compute the adjacency matrix.
            W = sparse( E(:,1), E(:,1), ones(p,1) );
            %symmetric
            W = max(W,W');
            %combinatorial Laplacian
            D = spdiags(sum(W)', 0, n,n);
            L = D-W;

        case 'contan'
        L = cotan_laplacian_matrix(vertex,faces) 
    end

end

