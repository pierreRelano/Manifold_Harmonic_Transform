path(path, 'mesh/');
path(path, 'mesh/meshes/');
path(path, 'mesh/toscahires-mat/');
path(path, 'mesh/mesh_deform/');
path(path, 'mesh/Models/');
path(path, 'lib/gptoolbox/mesh');
path(path, 'lib/gptoolbox/external/toolbox_fast_marching/toolbox');
path(path, 'lib/toolbox_graph/');% Mesh Smoothing  with Manifold Harmonics [Vallet et al. 2008]
% INPUTS:
% filePath: Path of the mesh you will use
% eigenNumber: Number of bases (Eigenvector of the Laplacian) you will use
% OUTPUTS:
% Displays original and modified version of the mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                               ^-_-^                                 %%%
%%%                    M. Cihan Ozer - UdeM, LIGUM                      %%%
%%%                        May 2016, Montreal 				%%%
%%%                                                                     %%%
%%% ALGORITHM:                                                          %%%
%%%                                                                     %%%
%%% 1. Assemble positive semi-definite discrete Laplacian               %%%
%%% 2. Compute eigenvectors for the Laplacian for getting MHB           %%%
%%% 3. Map the bases into canonical basis                               %%%
%%% 4. Transform the mesh into frequency space (MHT)                    %%%
%%% 5. Smooth the mesh                                                  %%%
%%% 6. Transform the mesh back into geometry space (MHT-1)              %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTES:
% - THIS SCRIPT DEPANDS ON GPTOOLBOX (Mesh reading/displaying and computation of 
%   cotangent Laplacian and Hodge star 0):
%   https://github.com/alecjacobson/gptoolbox
%
% - If you are using an old version of Matlab, ~ in operations (such as
%   [vertexNumber, ~] = size(vert)) may give you an error.
%   In this case you should replace ~ with a variable.
clc
clear all
close all


% READ MESH FILE
filePath = 'mesh/bunny_non_uniformly_sampled.off';
[~, remain] = strtok(filePath, '.');
if strcmp('.obj', remain) == 1
     [vertices, faces] = readOBJ(filePath);     %Read OBJ file
elseif strcmp('.off', remain) == 1
    [vertices, faces] = readOFF(filePath);  % Read OFF file
else
    disp('Unknown file type!');
end

% DISPLAY ORIGINAL MESH
figure('name', 'Before MHT');
plot_mesh(vertices, faces);
colormap gray(256)
shading faceted;
camlight;


% STEP 1: Assemble positive semi-definite discrete Laplacian
L = full(cotmatrix(vertices, faces));   % Get cotangent Laplacian
options.symmetrize = 1;
options.normalize = 0;
L_2 = compute_mesh_laplacian(vertex,face,'conformal',options);

pause
M = full(massmatrix(vertices, faces, 'barycentric'));   % Get dual area of vertices (Hodge star 0)
sprintf('Compupte Laplacian - done');   
Minv = sqrt(diag(1./diag(M)));    % Get M-1/2 for symmetry 
sprintf('Hodge O inversion - done'); 
% Minv1 = sqrt(inv(M));
%inv_diff = Minv - Minv1;
beltrami = Minv * L * Minv; % Get positive semi-definite discrete Laplacian (Eqauation 2)
beltrami = beltrami * -1;   % For positive eigenvalues
% Handle numerical precision issue:
% http://stackoverflow.com/a/33259074
beltrami = (beltrami + beltrami.') * 0.5;   % Now our Laplacian is symmetric, and
                                            % its eigenvectors are orthonormal
sprintf('Lapalce-Beltrami operator - done');                                            
% STEP 2: Compute eigenvectors for the Laplacian for getting MHB

% Apperantly, eig() function does not give orthogonal eigenvectors
% even if the Laplacian is positive semi-definite:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/29459 (First entry)
[~, eigen_val, eigen_vect] = svd(beltrami);
% nbr = 200;
% opts.disp = 0;
% [eigen_vect, eigen_val] = eigs(beltrami,nbr,'SM',opts);

% Sort eigenvectors by increasing eigenvalues (Ascending order)
[~, I] = sort(diag(eigen_val));
eigen_vect = eigen_vect(:, I);
sprintf('Eigen Decomposition - done'); 

% STEP 3: Map the eigen Vector bases into canonical basis
Hktemp = Minv * eigen_vect;
eigen_number = 250;
Hk = Hktemp(:,1:eigen_number);   % Take only the bases you will use


% Plot a sub-set of the eigenvectors on top of the topology
index_list = round(linspace(2,eigen_number,6));
tau=2.2; % saturation for display
clf;
figure('name', 'Subset of eigenvectors forming the the Manifold Harmonics Basis');
for i=1:length(index_list)
    v = real(Hk(:,index_list(i)));
    v = clamp(v/std(v),-tau,tau);
    options.face_vertex_color = v;
    subplot(2,3,i);
    plot_mesh(vertices,faces, options);
    shading interp; camlight; axis tight;
    colormap winter(256);
end

pause

% STEP 4: Transform the mesh into frequency space (MHT)
% For this operation, matrix multiplication is much faster than a loop
% This operation will give you a row vector
Xk = vertices(:,1)' * M * Hk;
Yk = vertices(:,2)' * M * Hk;
Zk = vertices(:,3)' * M * Hk;

% Construct MHT
MHT = zeros(eigen_number,3);
MHT(:,1) = Xk(1,:);
MHT(:,2) = Yk(1,:);
MHT(:,3) = Zk(1,:);

figure('name', 'MHB Spectrum');
plot(MHT); axis('tight');
legend('X', 'Y', 'Z');


% STEP 5: Design the filter
%smooth lowpass filter
[b_low,a_low] = butter(2,0.4,'low');
F_low = freqz(b_low,a_low,floor(eigen_number))';
hold on
plot((0:(1/(eigen_number - 1)):1), abs(F_low), 'r');

%sharp lowpass filter
[b_low2,a_low2] = butter(10,0.2,'low');
F_low2 = freqz(b_low2,a_low2,floor(eigen_number))';
hold on
plot((0:(1/(eigen_number - 1)):1), abs(F_low2), 'c');

%highpass
[b_high,a_high] = butter(2,0.2,'high');
F_high = freqz(b_high,a_high,floor(eigen_number))';
hold on
plot((0:(1/(eigen_number - 1)):1), abs(F_high), 'g');

%band reject filter 
[b_stop,a_stop] = butter(15,[0.3 0.7],'stop');
F_stop = freqz(b_stop,a_stop,floor(eigen_number))';
hold on
plot((0:(1/(eigen_number - 1)):1), abs(F_stop), 'm');

%design filter for band-exageration


% STEP 6: Smooth the mesh and transform it back into geometry space (MHT-1)
% Allocate memory for vertices
[vertexNumber, ~] = size(vertices); % Get vertex number
dumVertX = zeros(vertexNumber,1);
dumVertY = zeros(vertexNumber,1);
dumVertZ = zeros(vertexNumber,1);

%Choose filter from F_low, F_low2, F_high, F_stop
F = abs(F_low);

% MHT-1
for k = 1:eigen_number
     dumVertX = dumVertX + (F(1,k) * Xk(1,k) * Hk(:,k));
     dumVertY = dumVertY + (F(1,k) * Yk(1,k) * Hk(:,k));
     dumVertZ = dumVertZ + (F(1,k) * Zk(1,k) * Hk(:,k));
end

vertex_reconstruct = zeros(vertexNumber,3);
vertex_reconstruct(:,1) = dumVertX;
vertex_reconstruct(:,2) = dumVertY;
vertex_reconstruct(:,3) = dumVertZ;
% DISPLAY NEW MESH
figure('name', strcat('After MHT: ', int2str(eigen_number), ' bases is in use'));
plot_mesh(vertex_reconstruct, faces);
colormap gray(256)
shading faceted;
camlight;


% figure('name', 'Iterataive reconstructed meshes after the Inverse MHT adding Manifold Harmonics Basis components');
% index_list = round(linspace(2,eigenNumber,6));
% for i=1:length(index_list)
%     percent_modes_to_kept = index_list(i);
%     % set to zero high pass coefficients
%     q = round(percent_modes_to_kept/100*eigenNumber); % number of modes to filter
%     % set to 0 high-pass coeffs
%     Xk_temp = MHT(:,1)';
%     Yk_temp = MHT(:,2)';
%     Zk_temp = MHT(:,3)';
%     Xk_temp(q+1:end) = 0;
%     Yk_temp(q+1:end) = 0;
%     Zk_temp(q+1:end) = 0;
%     % MHT-1
%     for k = 1:eigenNumber
%          dumVertX = dumVertX + (Xk_temp(1,k) * Hk(:,k));
%          dumVertY = dumVertY + (Yk_temp(1,k) * Hk(:,k));
%          dumVertZ = dumVertZ + (Zk_temp(1,k) * Hk(:,k));
%     end
%     % MAP NEW VERTEX POSITIONS
%     Vfinal = zeros(vertexNumber,3);
%     Vfinal(:,1) = dumVertX;
%     Vfinal(:,2) = dumVertY;
%     Vfinal(:,3) = dumVertZ;
%     % display the reconstructed mesh
%     subplot(2,3,i);
%     plot_mesh(Vfinal,faces);
%     shading interp; camlight; axis tight;
%     colormap gray(256)
% end




