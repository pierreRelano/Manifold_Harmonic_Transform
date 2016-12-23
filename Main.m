path(path, 'meshes/');
path(path, 'lib/gptoolbox/mesh');
path(path, 'lib/gptoolbox/external/toolbox_fast_marching/toolbox');
path(path, 'lib/toolbox_graph/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ALGORITHM:                                                          %%%
%%% Mesh Filtering with Manifold Harmonics [Vallet et al. 2008]         %%%
%%%                                                                     %%%
%%%                                                                     %%%
%%% 1. Assemble positive semi-definite discrete Laplacian               %%%
%%% 2. Compute eigenvectors for the Laplacian for getting MHB           %%%
%%% 3. Map the bases into canonical basis                               %%%
%%% 4. Transform the mesh into frequency space (MHT)                    %%%
%%% 5. Smooth the mesh                                                  %%%
%%% 6. Transform the mesh back into geometry space (MHT-1)              %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;


% read mesh file
filePath = 'meshes/bunny_non_uniformly_sampled.off';
[~, remain] = strtok(filePath, '.');
if strcmp('.obj', remain) == 1
     [vertices, faces] = readOBJ(filePath);     %Read OBJ file
elseif strcmp('.off', remain) == 1
    [vertices, faces] = readOFF(filePath);  % Read OFF file
else
    disp('Unknown file type!');
end

%-------------------------------------------------------------
%               Step 0: Display original mesh
%-------------------------------------------------------------

figure('name', 'Before MHT');
plot_mesh(vertices, faces);
colormap gray(256)
%shading interp; axis tight;
shading faceted;
camlight;

%-------------------------------------------------------------
%               STEP 1: Build discrete Laplacian
%-------------------------------------------------------------


%L = full(cotmatrix(vertices, faces));   % Get cotangent Laplacian weight
%L2 = full(compute_cotangent_matrix(vertices, faces));
%L_diff = L - L2;
L = compute_manifold_laplacian(vertices, faces, 'combinatorial');
%options.symmetrize = 1;
%options.normalize = 0;
%L_2 = compute_mesh_laplacian(vertices,faces,'conformal',options);

msg = sprintf('Lapalce matrix - done');  
disp(msg);

M = full(massmatrix(vertices, faces, 'barycentric'));   % Get dual area of vertices (Hodge star 0)
msg = sprintf('Compupte Laplacian - done');   
disp(msg);
Minv = sqrt(diag(1./diag(M)));    % Get M-1/2 for symmetry / faster than inv()
msg = sprintf('Hodge O inversion - done'); 
disp(msg);
% Minv1 = sqrt(inv(M));
%inv_diff = Minv - Minv1;
beltrami = Minv * L * Minv; % Get positive semi-definite discrete Laplacian (Eqauation 2)
beltrami = beltrami * -1;   % For positive eigenvalues
% Handle numerical precision issue:
% http://stackoverflow.com/a/33259074
beltrami = (beltrami + beltrami.') * 0.5;   % Now our Laplacian is symmetric, and
                                            % its eigenvectors are orthonormal
msg = sprintf('Lapalce-Beltrami operator - done');  
disp(msg);

beltrami = (L + L.') * 0.5;

%-------------------------------------------------------------
%          STEP 2: Compute Eigenvectors Laplcian MHB
%-------------------------------------------------------------

% Apperantly, eig() function does not give orthogonal eigenvectors
% even if the Laplacian is positive semi-definite:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/29459 (First entry)
[~, eigen_val, eigen_vect] = svd(beltrami);
%nbr = 5000;
%opts.disp = 0;
%[eigen_vect, eigen_val] = eig(beltrami); %,nbr,'SM',opts);

% Sort eigenvectors by increasing eigenvalues (Ascending order)
[~, I] = sort(diag(eigen_val));
eigen_vect = eigen_vect(:, I);

msg = sprintf('Eigen Decomposition - done'); 
disp(msg);

%-------------------------------------------------------------
% STEP 3: Map the eigen Vector bases into canonical basis
%-------------------------------------------------------------

Hktemp = Minv * eigen_vect;
eigen_nmbr = 5000;
Hk = Hktemp(:,1:eigen_nmbr);   % Take only the bases you will use


% Plot a sub-set of the eigenvectors on top of the topology
%eigen_number_display = 250;
%index_list = round(linspace(2,eigen_number_display,6));
eigen_nmbr_list = [2, 50, 100, 250, 500, 1000];
tau=2.2; % saturation for display
clf;
figure('name', 'Subset of eigenvectors forming the the Manifold Harmonics Basis');
for i=1:length(eigen_nmbr_list)
    v = real(Hk(:,eigen_nmbr_list(i)));
    v = clamp(v/std(v),-tau,tau); %clamp a value- function defined in gptoolbox
    options.face_vertex_color = v;
    subplot(2,3,i);
    plot_mesh(vertices,faces, options);
    shading interp; camlight; axis tight;
    colormap winter(256);
    str=sprintf('m = %d', eigen_nmbr_list(i));
    title(str);
end

msg = sprintf('MHB - done'); 
disp(msg);

%-------------------------------------------------------------
%    STEP 4: Transform the mesh into frequency space (MHT)
%-------------------------------------------------------------

% For the projection, matrix multiplication is much faster than a loop
% This operation will give you three row vectors
Xk = vertices(:,1)' * M * Hk;
Yk = vertices(:,2)' * M * Hk;
Zk = vertices(:,3)' * M * Hk;

% Construct MHT
MHT = zeros(eigen_nmbr,3);
MHT(:,1) = Xk(1,:);
MHT(:,2) = Yk(1,:);
MHT(:,3) = Zk(1,:);

figure('name', 'MHB Spectrum');
plot(MHT); axis('tight');
legend('X', 'Y', 'Z');

msg = sprintf('MHT - done'); 
disp(msg);

%-------------------------------------------------------------
%               STEP 5: Design the filters
%-------------------------------------------------------------

figure('name', 'Frequency filters');
%smooth lowpass filter
[b_low,a_low] = butter(2,0.2,'low');
F_low = freqz(b_low,a_low,floor(eigen_nmbr))';
hold on
plot((0:(1/(eigen_nmbr - 1)):1), abs(F_low), 'r');

%sharp lowpass filter
[b_low2,a_low2] = butter(10,0.2,'low');
F_low2 = freqz(b_low2,a_low2,floor(eigen_nmbr))';
hold on
plot((0:(1/(eigen_nmbr - 1)):1), abs(F_low2), 'c');

%highpass
[b_high,a_high] = butter(2,0.2,'high');
F_high = freqz(b_high,a_high,floor(eigen_nmbr))';
hold on
plot((0:(1/(eigen_nmbr - 1)):1), abs(F_high), 'g');

%high freq exageration
[b_high_exa,a_high_exa] = butter(4,0.95,'high');
F_high_exa = freqz(b_high_exa,a_high_exa,floor(eigen_nmbr))';
F_high_exa_abs = abs(F_high_exa)*7+1;
hold on
plot((0:(1/(eigen_nmbr - 1)):1), F_high_exa_abs, '--b');

%low freq exageration
[b_low_exa,a_low_exa] = butter(6,0.5,'low');
F_low_exa = freqz(b_low_exa,a_low_exa,floor(eigen_nmbr))';
F_low_exa_abs = abs(F_low_exa)*10 +1; %FAT exageration!
hold on
plot((0:(1/(eigen_nmbr - 1)):1), F_low_exa_abs, 'c');

%band reject filter 
[b_stop,a_stop] = butter(15,[0.3 0.7],'stop');
F_stop = freqz(b_stop,a_stop,floor(eigen_nmbr))';
hold on
plot((0:(1/(eigen_nmbr - 1)):1), abs(F_stop), 'm');

%design bandpass filter
[b_stop,a_stop] = butter(3,[0.3 0.7],'bandpass');
F_stop = freqz(b_stop,a_stop,floor(eigen_nmbr))';
hold on
plot((0:(1/(eigen_nmbr - 1)):1), abs(F_stop), 'b');

legend('smooth lowpass', 'sharp lowpass', 'highpass', ...
    'high freqs exageration', 'low freqs exageration', ...
    'band-reject', 'bandpass');

%-------------------------------------------------------------
%   STEP 6: Filter mesh and inverse MHT into geometry space 
%-------------------------------------------------------------

% Allocate memory for vertices
[vertex_nmbr, ~] = size(vertices); % Get vertex number
IMHT_vert_X = zeros(vertex_nmbr,1);
IMHT_vert_Y = zeros(vertex_nmbr,1);
IMHT_vert_Z = zeros(vertex_nmbr,1);


%=========================
% Perfect reconstruction
%=========================

for k = 1:eigen_nmbr
     IMHT_vert_X = IMHT_vert_X + (Xk(1,k) * Hk(:,k));
     IMHT_vert_Y = IMHT_vert_Y + (Yk(1,k) * Hk(:,k));
     IMHT_vert_Z = IMHT_vert_Z + (Zk(1,k) * Hk(:,k));
end
% assign new vertex position
vertex_reconstruct = zeros(vertex_nmbr,3);
vertex_reconstruct(:,1) = IMHT_vert_X;
vertex_reconstruct(:,2) = IMHT_vert_Y;
vertex_reconstruct(:,3) = IMHT_vert_Z;


%=========================
% Filtered reconstruction
%=========================
%
% %Choose filter from F_low, F_low2, F_high, F_stop
% F = F_low_exa_abs; %abs(F_low);
% 
% % MHT-1
% for k = 1:eigen_number
%      IMHT_vert_X = IMHT_vert_X + (F(1,k) * Xk(1,k) * Hk(:,k));
%      IMHT_vert_Y = IMHT_vert_Y + (F(1,k) * Yk(1,k) * Hk(:,k));
%      IMHT_vert_Z = IMHT_vert_Z + (F(1,k) * Zk(1,k) * Hk(:,k));
% end
% % assign new vertex position
% vertex_reconstruct = zeros(vertexNumber,3);
% vertex_reconstruct(:,1) = IMHT_vert_X;
% vertex_reconstruct(:,2) = IMHT_vert_Y;
% vertex_reconstruct(:,3) = IMHT_vert_Z;


%=========================
% Iterative reconstruction
%=========================
% 
% figure('name', 'Iterataive reconstructed meshes after the Inverse MHT adding Manifold Harmonics Basis components');
% eigen_nmbr_list = [6, 10, 25, 50, 100, 250, 500, 1000, 5000];
% for i=1:length(eigen_nmbr_list)
%     eigen_number_kept = eigen_nmbr_list(i);
%     IMHT_vert_X = zeros(vertex_nmbr,1);
%     IMHT_vert_Y = zeros(vertex_nmbr,1);
%     IMHT_vert_Z = zeros(vertex_nmbr,1);
%     % MHT-1
%     for k = 1:eigen_number_kept
%         IMHT_vert_X = IMHT_vert_X + (Xk_temp(1,k) * Hk(:,k));
%         IMHT_vert_Y = IMHT_vert_Y + (Yk_temp(1,k) * Hk(:,k));
%         IMHT_vert_Z = IMHT_vert_Z + (Zk_temp(1,k) * Hk(:,k));
%     end
%     % assign new vertex position
%     vertex_reconstruct = zeros(vertex_nmbr,3);
%     vertex_reconstruct(:,1) = IMHT_vert_X;
%     vertex_reconstruct(:,2) = IMHT_vert_Y;
%     vertex_reconstruct(:,3) = IMHT_vert_Z;
%     % display the reconstructed mesh
%     subplot(3,3,i);
%     plot_mesh(vertex_reconstruct,faces);
%     shading interp; camlight; axis tight;
%     colormap gray(256)
%     str=sprintf('m = %d', eigen_number_kept);
%     title(str);
% end

msg = sprintf('IMHT - done'); 
disp(msg);

%=============================================
% Display reconstructed mesh vs original mesh
%=============================================

figure('name', 'Original and reconstructed mesh');
subplot(1,2,1);
plot_mesh(vertices, faces);
title('Orignal Mesh')
shading interp; axis tight;
subplot(1,2,2);
plot_mesh(vertex_reconstruct, faces);
title(strcat('Reconstruction after MHT-IMHT (m = 5000)'));
colormap gray(256)
shading interp; axis tight;


%save workspace_mesh_independent_laplacian.mat
