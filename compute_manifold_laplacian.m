function L = compute_manifold_laplacian(vertices, faces, type)
%The combinatorial laplacian is a linear operator (thus a NxN matrix where N is the number of vertices). 
%It depends only on the connectivity of the mesh, thus on face only.

    %TODO: test the indice of combinatorial laplacian
    switch lower(type)
        case 'combinatorial'
            n = length(vertices(:,1));
            %Compute Adjacency matrix
            W = adjacency_matrix(faces);
            %Degree matrix
            D = full(spdiags(sum(W)', 0, n,n));
            %symmetric combinatorial Laplacian
            L = D-W;
        case 'contan'
        L = cotan_laplacian_matrix(vertices,faces); 
    end

end

