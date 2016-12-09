path(path, 'lib/toolbox_graph/');
path(path, 'lib/toolbox_graph/off/');
path(path, 'lib/toolbox_graph/toolbox/');
path(path, 'mesh/meshes/');
path(path, 'mesh/toscahires-mat/');
path(path, 'mesh/mesh_deform/');
path(path, 'lib/toolbox_signal/');
path(path, 'lib/toolbox_general/');

clear
clc
clear options;

%------------------------
%load mesh

name = 'elephant-50kv';
[vertex,faces] = read_mesh([name '.off']);
n = size(vertex,2);
options.name = name;


%------------------------
%display cutting limit mesh

W = ones(n,1);
W(vertex(1,:)<mean(vertex(1,:))) = .4;
options.W = W;
clf;
hold on;
options.face_vertex_color = W;
plot_mesh(vertex,faces, options);
colormap jet(256);

%------------------------
%Perform front propagation using this speed function.

landmarks = [500];
options.constraint_map = [];
[D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks, options);

%Display the distance map.
clf;
hold on;
options.face_vertex_color = mod( 20*D/max(D),1 );
plot_mesh(vertex,faces, options);
colormap jet(256);
h = plot3(vertex(1,landmarks), vertex(2,landmarks), vertex(3,landmarks), 'r.');
set(h, 'MarkerSize', 20);


%------------------------
%downsampled mesh

% uniformly sampled
options.W = ones(n,1);
%%non-uniformly sampled
% W = ones(n,1);
% W(vertex(1,:)<mean(vertex(1,:))) = .4;
% options.W = W;

m = 2500;
clf;
k = 1; displist = (2500);
for i=2:m
    % select
    [tmp,landmarks(end+1)] = max(D);
    % update
    options.constraint_map = D;
    [D1,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks,options);
    D = min(D,D1);
    if i==displist(k)
        i
        [D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks);
        % compute the mesh
        V = Q(faces); V = sort(V,1);
        V = unique(V', 'rows')';
        d = 1 + (V(1,:)~=V(2,:)) + (V(2,:)~=V(3,:));
        %
        I = find(d==3); I = sort(I);
        z = zeros(n,1);
        z(landmarks) = (1:length(landmarks))';
        facesV = z(V(:,I));
        vertexV = vertex(:,landmarks);
        % Re-orient the faces so that they point outward of the mesh.
        options.method = 'slow';
        options.verb = 0;
        facesV = perform_faces_reorientation(vertexV,facesV);
        write_off([name '_uniformly_sampled.off'],vertexV, facesV)
        % display
        subplot(1,2,k);
        options.face_vertex_color = [];
        plot_mesh(vertexV,facesV, options);
        shading faceted;
        %
        k = k+1;
    end
end

pause

%------------------------
% bi-partite noised mesh

%Right part noise
noise = randn(size(vertex))*.002;
noise(:,vertex(1,:)<mean(vertex(1,:))) = 0;
vertex1 = vertex+noise;
write_off('bunny_half_noise.off',vertex1, faces)

clf;
figure('name', 'Before MHT');
plot_mesh(vertex1,faces);
colormap gray(256)
% add the display of the triangle
%shading faceted;
% when you rotate the camera, focuss the light using
camlight;

%------------------------
% Add noise to the full mesh

% compute the normals
normals = compute_normal(vertex,faces);
rho = randn(1,size(vertex,2))*.002;
vertex_noise = vertex + repmat(rho,3,1).*normals;

% display the mesh
clf;
figure('name', 'Before MHT');
plot_mesh(vertex_noise,faces);
colormap gray(256)
% add the display of the triangle
%shading faceted;
% when you rotate the camera, focuss the light using
camlight;
