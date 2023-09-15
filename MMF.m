clear all

P = BPMmatlab.model;

% This example shows that modes can be found and stored in the parameters
% struct using the findModes function, then used as initial E fields. Here,
% we use the LP31o-like mode of a non-trivial refractive index profile.

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs

% %% Visualization parameters
% P.updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state.

%% Resolution-related parameters (check for convergence)
P.Lx_main = 200e-6;        % [m] x side length of main area
P.Ly_main = 200e-6;        % [m] y side length of main area
P.Nx_main = 200;          % x resolution of main area
P.Ny_main = 200;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 5e-7; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1550e-9; % [m] Wavelength
NA = 0.2;
n_clad = 1.444; % pure silica at 1550nm
n_exp = sqrt(NA^2 + n_clad^2);
V = (2*pi*25e-6*NA)/P.lambda; % V^2/8 = 51 modes in graded index
P.n_background = 1.444; % MMF do not have capillary, assumingm all cladding outside of the MMF core
n_0 = sqrt(n_clad^2 + (V*P.lambda/(2*pi*25e-6))^2); % fibre core,(Ge doped silica)  assume V = 20.2683,n_0 = n_core = 1.4578
P.n_0 = n_0; %n_core = 1.4578


[X,Y] = ndgrid(single(P.x),single(P.y));
P.n.n = single(calcRI(X,Y,n_clad,P.n_background,n_0));
P.n.Lx = P.dx*size(P.n.n,1);
P.n.Ly = P.dy*size(P.n.n,2);
P.n.xSymmetry = P.xSymmetry;
P.n.ySymmetry = P.ySymmetry;

% We search for the 10 modes with effective refractive index closest to
% n_0:
P = findModes(P,7);
save('MMF_n143.mat', 'P');

% 
% for k1 = 1:10
%     for k2 = 1:10
%         oi(k1,k2) = sum(sum(conj(SMF.P.modes(k1).field).*P.modes(k2).field));  %conjugate SMF and MMF 
%     end
%     
% end
% figure()
% imagesc(abs(oi))




figure()
mesh(P.n.n)
ylabel('Y-Position(μm)');
xlabel('X-Position(μm)');
title('Cross section of the MMF')  %

figure()
mesh(real(P.modes(1).field) );view([ 0 0 1])

figure()
subplot(1,2,1)
mesh(real(P.modes(2).field) );view([ 0 0 1])
subplot(1,2,2)
mesh(real(P.modes(3).field) );view([ 0 0 1])

% P.modes.neff

% % Run solver
% FD_BPM(P);

%% USER DEFINED RI FUNCTIONS
function n = calcRI(X,Y,n_cladding,n_background,n_core)
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background

% Cladding 1 is step index(center):
n((X + 0e-6).^2 + (Y + 0e-6).^2 < 62.5e-6^2) = n_cladding;

% Core 1 is graded index:
corepos = [0 0];
r = 25e-6;
R = sqrt((X - corepos(1)).^2 + (Y - corepos(2)).^2);
n(R < r) = n_background + (n_core - n_background)*(1 - (R(R < r)/r).^2); % Equation for parabolic graded index core

end
