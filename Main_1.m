%% STEP 0: Set workspace path
path(path, 'lib/toolbox_graph/');
path(path, 'lib/toolbox_graph/off/');
path(path, 'lib/toolbox_graph/toolbox/');
path(path, 'mesh/');
path(path, 'mesh/meshes/');
path(path, 'mesh/toscahires-mat/');
path(path, 'mesh/mesh_deform/');
path(path, 'lib/gptoolbox/mesh/');

%% ---------------------------------------------------------
%            STEP 1: Load and display a Manifold
% ----------------------------------------------------------

% large meshes, you can try 'skull', 'fandisk', 'lion-head', 'bunny'
% small meshes, you can try 'mushroom', 'nefertiti', 'hammerheadtriang'
clear
clc

name_list = {'bunny'; 'bunny_non_uniformly_sampled'; 'bunny_noise'; 'null0'; 'null0_uniformly_sampled'; 'venus'};
name = 'venus';

% load from file
[vertex,face] = readOFF([name '.off']);
nface = size(face,1);
%nvert = max(max(face));

% display the mesh
clf;
figure('name', 'Before MHT');
plot_mesh(vertex,face);
%view(-[0 0]);
colormap gray(256)
% add the display of the triangle
shading faceted;
% when you rotate the camera, focuss the light using
camlight;


%% ---------------------------------------------------------
%            STEP 2: Compute discrete Laplacians 
% ----------------------------------------------------------

%combinatorial Laplcian computation
%L = compute_combinatorial_manifold_laplacian(face);
% 
%%cotangent Laplcian computation
% options.symmetrize = 1;
% options.normalize = 0;
% L = compute_mesh_laplacian(vertex,face,'conformal',options);
%L_noise = compute_mesh_laplacian(vertex_noise,face,'conformal',options);

L = full(cotmatrix(vertex,face));   % Get cotangent Laplacian
M = full(hodge_star_0(vertex, face));   % Get dual area of vertices (Hodge star 0)
Minv = sqrt(inv(M)); % Get M-1/2 for symmetry 
beltrami = Minv * L * Minv; % Get positive semi-definite discrete Laplacian (Eqauation 2)
beltrami = beltrami * -1;   % For positive eigenvalues
% Handle numerical precision issue:
% http://stackoverflow.com/a/33259074
beltrami = (beltrami + beltrami.') * 0.5;   % Now our Laplacian is symmetric, and
                                            % its eigenvectors are orthonormal
sprintf('Compupte Laplacian - done')   

% STEP 2: Compute eigenvectors for the Laplacian for getting MHB

% Apperantly, eig() function does not give orthogonal eigenvectors
% even if the Laplacian is positive semi-definite:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/29459 (First entry)
[~, e, eVec] = svd(beltrami);

% Sort eigenvectors by increasing eigenvalues (Ascending order)
[~, I] = sort(diag(e));
eVec = eVec(:, I);

sprintf('Eigenstuff - done')

% STEP 3: Map the bases into canonical basis
Hktemp = Minv * eVec;
eigenNumber = size(Hktemp,1);
Hk = Hktemp(:,1:eigenNumber);   % Take only the bases you will use

% STEP 4: Transform the mesh into frequency space (MHT)
% For this operation, matrix multiplication is much faster than a loop
% This operation will give you a row vector
Xk = vertex(:,1)' * M * Hk;
Yk = vertex(:,2)' * M * Hk;
Zk = vertex(:,3)' * M * Hk;

% STEP 5-6: Smooth the mesh and transform it back into geometry space (MHT-1)


% amplify some band pass coefficients
q = 10; % number of modes to filter
% set to 0 high-pass coeffs
Xk(1, q+1:q+q) = Xk(1, q+1:q+q).*.1;
Yk(1, q+1:q+q) = Yk(1, q+1:q+q).*1.1;
Zk(1, q+1:q+q) = Zk(1, q+1:q+q).*1.1;


% Allocate memory for vertices
[vertexNumber, ~] = size(vertex); % Get vertex number
dumVertX = zeros(vertexNumber,1);
dumVertY = zeros(vertexNumber,1);
dumVertZ = zeros(vertexNumber,1);

% MHT-1
for k = 1:eigenNumber
     dumVertX = dumVertX + (Xk(1,k) * Hk(:,k));
     dumVertY = dumVertY + (Yk(1,k) * Hk(:,k));
     dumVertZ = dumVertZ + (Zk(1,k) * Hk(:,k));
end

% MAP NEW VERTEX POSITIONS
Vfinal = zeros(vertexNumber,3);
Vfinal(:,1) = dumVertX;
Vfinal(:,2) = dumVertY;
Vfinal(:,3) = dumVertZ;

% DISPLAY NEW MESH
figure('name', strcat('After MHT: ', int2str(eigenNumber), ' bases is in use'));
plot_mesh(Vfinal, face);
shading interp; camlight; axis tight;
colormap gray(256)

stop

%% ---------------------------------------------------------
%           STEP 3: Performing eigendecomposition
% ----------------------------------------------------------


if size(vertex,2)<1000
    [MHB,lambda] = eig(full(L));
     nbr = size(MHB,1);
else
    nbr = 200;
    opts.disp = 0;
    [MHB,lambda] = eigs(L,nbr,'SM',opts);
end
% order the eigenvector by increasing frequencies
lambda = diag(lambda);
[lambda_ordered,Index] = sort(lambda, 'ascend');
MHB = real(MHB(:,Index));

% plot the eigenvalues, spectrum of the triangulation. 
% depends only on the topology of the mesh
clf;
figure('name', 'Eigenvalues, spectrum of the triangulation');
plot(lambda_ordered); axis('tight');

% plot a sub-set of the eigenvectors on top of the topology
index_list = round(linspace(2,nbr, 6));
tau=2.2; % saturation for display
clf;
figure('name', 'Subset of eigenvectors forming the the Manifold Harmonics Basis');
for i=1:length(index_list)
    v = real(MHB(:,index_list(i)));
    v = clamp(v/std(v),-tau,tau);
    options.face_vertex_color = v;
    subplot(2,3,i);
    plot_mesh(vertex,face);
    shading interp; camlight; axis tight;
    colormap winter(256);
end

% % extract one of the eigenvectors
% c = MHB(:,1);
% % assign a color to each vertex
% options.face_vertex_color = rescale(c, 0,255);
% % display
% clf;
% plot_mesh(vertex,face, options);

sprintf('Eigendecomposition - done')  
%% ---------------------------------------------------------
%                    STEP 4: Design filter
% ----------------------------------------------------------



%% ---------------------------------------------------------
%                  STEP 5: MHT and Filtering
% ----------------------------------------------------------

% projection of the geometry in the laplacian basis
MHT = (MHB' * vertex'); 

%plot the spectrum
clf;
plot(MHT); axis('tight');
legend('X', 'Y', 'Z');
figure('name', 'Iterataive reconstructed meshes after the Inverse MHT with more and more Manifold Harmonics components');
index_list = round(linspace(2,100,6));

percent_modes_to_kept = index_list(i);
% amplify some band pass coefficients
q = 10; % number of modes to filter
% set to 0 high-pass coeffs
MHT(q+1:q+q,:) = MHT(q+1:q+q,:).*0.1;
% Inverse MH
vertex1 = (MHB * MHT)';
% display the reconstructed mesh
plot_mesh(vertex1,face);
shading interp; camlight; axis tight;
colormap gray(256)


% for i=1:length(index_list)
%     if i ~= 1 
%         MHT = (MHB' * vertex');
%     end
%     percent_modes_to_kept = index_list(i);
%     % set to zero high pass coefficients
%     q = round(percent_modes_to_kept/100*nbr); % number of modes to filter
%     % set to 0 high-pass coeffs
%     MHT(q+1:end,:) = 0;
%     % Inverse MHT
%     vertex1 = (MHB * MHT)';
%     % display the reconstructed mesh
%     subplot(2,3,i);
%     plot_mesh(vertex1,face);
%     shading interp; camlight; axis tight;
%     colormap gray(256)
% 
% end
