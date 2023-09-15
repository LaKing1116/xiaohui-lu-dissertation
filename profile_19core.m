clear all
M1 = zeros(1,10);
P = BPMmatlab.model;

% This code is for the refractive index profile of 19 SMFs

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs

% %% Visualization parameters
% P.updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state.

%% Resolution-related parameters (check for convergence)
P.Lx_main = 400e-6;        % [m] x side length of main area  
P.Ly_main = 400e-6;        % [m] y side length of main area
P.Nx_main = 400;          % x resolution of main area
P.Ny_main = 400;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 5e-7; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1550e-9; % [m] Wavelength
n_clad = 1.444; % pure silica at 1550nm
P.n_background = 1.42; % [] (may be complex) Background refractive index,  change n_background in this file and MMF_out file 
n_0 = sqrt(n_clad^2 + (2.3*P.lambda/(2*pi*2.5e-6))^2); % [] reference refractive index (fibre core) if assume V = 2.3  n_0=n_core=1.4617
P.n_0 = n_0;                                           
capillary_delta_before_tapering = (n_0^2-n_clad^2)/(2*(n_0^2))% 0.0121
capillary_delta_after_tapering = (n_clad^2-P.n_background^2)/(2*(n_clad^2))% 0.0096483, n_clad becomes new core, n_background becomes new cladding

[X,Y] = ndgrid(single(P.x),single(P.y));
P.n.n = single(calcRI(X,Y,n_clad,P.n_background,n_0));
P.n.Lx = P.dx*size(P.n.n,1);
P.n.Ly = P.dy*size(P.n.n,2);
P.n.xSymmetry = P.xSymmetry;
P.n.ySymmetry = P.ySymmetry;

% We search for the 10 modes with effective refractive index closest to
% n_0:
% P = findModes(P,10);

figure()
mesh(P.n.n)
ylabel('Vertical spatial distance(μm)');
xlabel('Horizontal spatial distance(μm)');
title('Cross section of nineteen SMFs')  %

figure()
mesh(real(P.modes(1).field) );view([ 0 0 1])

figure()
subplot(1,2,1)
mesh(real(P.modes(2).field) );view([ 0 0 1])
subplot(1,2,2)
mesh(real(P.modes(3).field) );view([ 0 0 1])


%% USER DEFINED RI FUNCTIONS
function n = calcRI(X,Y,n_cladding,n_background,n_core)
% n may be complex ,r = 80 / (2*sin(pi/12));
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
% cladding 1 is step index:(center)
n((X + 0e-6).^2 + (Y + 0e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2 (bottom)
n((X + 80e-6).^2 + (Y + 0e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 3 (top)
n((X - 80e-6).^2 + (Y + 0e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 4 (bottom-left)
n((X + 40e-6).^2 + (Y + 3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 5 (top-left)
n((X - 40e-6).^2 + (Y + 3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 6 (bottom-right)
n((X + 40e-6).^2 + (Y - 3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 7 (top-right)
n((X - 40e-6).^2 + (Y - 3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2-1 ()
n((X + 2*80e-6).^2 + (Y + 0e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2-2 ()
n((X - 2*80e-6).^2 + (Y + 0e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2-3 ()
n((X + 2*40e-6).^2 + (Y + 2*3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2-4 ()
n((X - 2*40e-6).^2 + (Y + 2*3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2-5 ()
n((X + 2*40e-6).^2 + (Y - 2*3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2-6 ()
n((X - 2*40e-6).^2 + (Y - 2*3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2-7 ()
n((X + 0e-6).^2 + (Y + 2*3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2-8 ()
n((X + 0e-6).^2 + (Y - 2*3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2-9 ()
n((X + 3^0.5*3^0.5*40e-6).^2 + (Y + 3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2-10 ()
n((X + 3^0.5*3^0.5*40e-6).^2 + (Y - 3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2-11 ()
n((X - 3^0.5*3^0.5*40e-6).^2 + (Y - 3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2-12 ()
n((X - 3^0.5*3^0.5*40e-6).^2 + (Y + 3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;






% Core 1 is step index:(center)
n((X + 0e-6).^2 + (Y + 0e-6).^2 < 2.5e-6^2) = n_core;

% % Core 2 (bottom)
n((X + 80e-6).^2 + (Y + 0e-6).^2 < 2.5e-6^2) = n_core;

% % Core 3 (top)
n((X - 80e-6).^2 + (Y + 0e-6).^2 < 2.5e-6^2) = n_core;

% % Core 4 (bottom-left)
n((X + 40e-6).^2 + (Y + 3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % Core 5 (top-left)
n((X - 40e-6).^2 + (Y + 3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % Core 6 (bottom-right)
n((X + 40e-6).^2 + (Y - 3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % Core 7 (top-right)
n((X - 40e-6).^2 + (Y - 3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;


% % core 2-1 ()
n((X + 2*80e-6).^2 + (Y + 0e-6).^2 < 2.5e-6^2) = n_core;

% % core 2-2 ()
n((X - 2*80e-6).^2 + (Y + 0e-6).^2 < 2.5e-6^2) = n_core;

% % core 2-3 ()
n((X + 2*40e-6).^2 + (Y + 2*3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % core 2-4 ()
n((X - 2*40e-6).^2 + (Y + 2*3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % core 2-5 ()
n((X + 2*40e-6).^2 + (Y - 2*3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % core 2-6 ()
n((X - 2*40e-6).^2 + (Y - 2*3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % core 2-7 ()
n((X + 0e-6).^2 + (Y + 2*3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % core 2-8 ()
n((X + 0e-6).^2 + (Y - 2*3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % core 2-9 ()
n((X + 3^0.5*3^0.5*40e-6).^2 + (Y + 3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % core 2-10 ()
n((X + 3^0.5*3^0.5*40e-6).^2 + (Y - 3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % core 2-11 ()
n((X - 3^0.5*3^0.5*40e-6).^2 + (Y - 3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % core 2-12 ()
n((X - 3^0.5*3^0.5*40e-6).^2 + (Y + 3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

end
