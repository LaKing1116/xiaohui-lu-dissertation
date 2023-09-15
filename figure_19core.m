clear all

E_normalized_SMF = zeros(600,600);
E_normalized_MMF = zeros(600,600);
P = BPMmatlab.model;

% This code is for evolution of modes throughout the tapered transition of the photonic lantern with 19 SMFs.

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
NA = 0.2;
n_clad = 1.444; % pure silica at 1550nm
n_exp = sqrt(NA^2 + n_clad^2);
V = (2*pi*2.5e-6*NA)/P.lambda;
capillary_delta = -0.01;
% P.n_background = 1.43; % Background refractive index  F doped silica,capillary
P.n_background = sqrt((-n_clad^2)/(2*capillary_delta-1)); % n_Background refractive index = 1.4298;F doped silica,capillary
n_0 = sqrt(n_clad^2 + (V*P.lambda/(2*pi*2.5e-6))^2); % [] reference refractive index (fibre core)(Ge doped silica)  assume V = 2.0268,n_0=n_core=1.4578
P.n_0 = n_0;                                         
capillary_delta_before_tapering = (n_0^2-n_clad^2)/(2*(n_0^2));% 0.0094
capillary_delta_after_tapering = (P.n_background^2-n_clad^2)/(2*(P.n_background^2));% -0.0098, n_clad becomes new core, n_background becomes new cladding

P.n.xSymmetry = P.xSymmetry;
P.n.ySymmetry = P.ySymmetry;

taper_ratio = [2 5 8 11 14 17];  %[2 5 8 11 14 17]
for k1 = 1:length(taper_ratio)

[X,Y] = ndgrid(single(P.x),single(P.y));
P.n.n = single(calcRI(X,Y,n_clad,P.n_background,n_0));
P.n.Lx = P.dx*size(P.n.n,1)/taper_ratio(k1); %%The taper is being achieved by re-scaling the physical coordinates of your simulations.
P.n.Ly = P.dy*size(P.n.n,2)/taper_ratio(k1); %%The taper ratio is the ratio between the diameter of the capillary before tapering and the diameter of the capillary after tapering.

P = findModes(P,51);  %MMF has 51 modes ; MMF's V^2/8 = 51.35

MMF=load('MMF_19.mat');

neff(k1,:) = [P.modes.neff];    %find neff after taper

for k3 = 1:19
    for k2 = 1:19
        oi(k3,k2) = sum(sum(conj(P.modes(k3).field).*MMF.P.modes(k2).field));   %  A=[n x n]; sum (A)=sum for n columns; sum(sum(A)) = sum for n rows
        
        E_normalized_SMF = (P.modes(k3).field)/(sqrt(sum(sum(abs(P.modes(k3).field.^2))))); %normalise the mode field distribution for unitary power
        E_normalized_MMF = (MMF.P.modes(k2).field)/(sqrt(sum(sum(abs(MMF.P.modes(k2).field.^2))))); %normalization
        E_checknormalized_SMF(k3,k2) = abs(sum(sum(E_normalized_SMF.*E_normalized_SMF)));  %result = 1, a check to make sure the normalization was done correctly. 
        E_checknormalized_MMF(k3,k2) = abs(sum(sum(E_normalized_MMF.*E_normalized_MMF)));  %result = 1, a check to make sure the normalization was done correctly. 
    
    end
    
end

% figure()
% imagesc(abs(oi))
A6(1,k1)= sum(sum(abs(oi))); %Get the magnitude (or absolute value) of a complex number  ,coupling efficiency = sum of power terms in coupling integral matrix
% A(2,k1)= abs(trace(oi));
end

% Divide original outer diameter 5*80(μm) by each element in a taper_ratio matrix
B = 400 ./ taper_ratio;
figure()
plot(B,neff)  % 1./taper_ratio(outer diameter) ---x axis,  neff----y axis
xlabel('Outer Diameter(μm)')  % Sets the label for the x-axis
ylabel('Neff')  % Sets the label for the y-axis
title('Evolution of modes for PL with 19 SMFs and -1% capillary delta')  % Sets the title for the entire plot
% % Set the limits of the x and y axes
% xlim([20 80])  % Replace xmin and xmax with your actual limits
% ylim([1.435 1.445])  % Replace ymin and ymax with your actual limits

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
