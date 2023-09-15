clear all

P = BPMmatlab.model;

% This code is for ideal coupling matrix

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

taper_ratio = [6:0.1:6.1];
for k1 = 1:length(taper_ratio)

[X,Y] = ndgrid(single(P.x),single(P.y));
P.n.n = single(calcRI(X,Y,n_clad,P.n_background,n_0));
P.n.Lx = P.dx*size(P.n.n,1)/taper_ratio(k1); %taper transition
P.n.Ly = P.dy*size(P.n.n,2)/taper_ratio(k1);

P = findModes(P,7);

MMF=load('MMF_n143.mat');

neff(k1,:) = [P.modes.neff];    %find neff after taper

for k3 = 1:7
    for k2 = 1:7
        oi(k3,k2) = sum(sum(conj(P.modes(k3).field).*P.modes(k2).field));   %  A=[n x n]; sum (A)=sum for n columns; sum(sum(A)) = sum for n rows
        
        E_normalized_SMF = (P.modes(k3).field)/(sqrt(sum(sum(abs(P.modes(k3).field.^2))))); %normalise the mode field distribution for unitary power;normalization = E-SMF/sqrt(power)
        E_normalized_MMF = (MMF.P.modes(k2).field)/(sqrt(sum(sum(abs(MMF.P.modes(k2).field.^2))))); %normalization = E-MMF/sqrt(power)
        E_checknormalized_SMF(k3,k2) = abs(sum(sum(E_normalized_SMF.*E_normalized_SMF)));  %result = 1, a check to make sure the normalization was done correctly. 
        E_checknormalized_MMF(k3,k2) = abs(sum(sum(E_normalized_MMF.*E_normalized_MMF)));  %result = 1, a check to make sure the normalization was done correctly. 
    end
    
end
figure()
imagesc(abs(oi))
colorbar;
xlabel('MMF Mode');
ylabel('SMF Mode');
title('Coupling matrix for PL with 7 SMFs');
% Coupling Efficiency
A3(1,k1)= 100*((sum(sum(abs(oi).^2)))/(size(oi, 1))); %size(oi, 1) is total input power

%% CROSSTALK  teacher method  row(SMF)

oi = double(oi);
power  = abs(oi).^2; %E = A^2
xt_row = (sum(power,2)-max(power,[],2))./max(power,[],2);
xt_dB = 10*log10(xt_row);

end


%% USER DEFINED RI FUNCTIONS
function n = calcRI(X,Y,n_cladding,n_background,n_core)
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
% cladding 1 is step index:(center)
n((X + 0e-6).^2 + (Y + 0e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 2 (right)
n((X + 80e-6).^2 + (Y + 0e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 3 (left)
n((X - 80e-6).^2 + (Y + 0e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 4 (bottom-right)
n((X + 40e-6).^2 + (Y + 3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 5 (bottom-left)
n((X - 40e-6).^2 + (Y + 3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 6 (top-right)
n((X + 40e-6).^2 + (Y - 3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;

% % cladding 7 (top-left)
n((X - 40e-6).^2 + (Y - 3^0.5*40e-6).^2 < 40e-6^2) = n_cladding;


% Core 1 is step index:(center)
n((X + 0e-6).^2 + (Y + 0e-6).^2 < 2.5e-6^2) = n_core;

% % Core 2 (right)
n((X + 80e-6).^2 + (Y + 0e-6).^2 < 2.5e-6^2) = n_core;

% % Core 3 (left)
n((X - 80e-6).^2 + (Y + 0e-6).^2 < 2.5e-6^2) = n_core;

% % Core 4 (bottom-right)
n((X + 40e-6).^2 + (Y + 3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % Core 5 (bottom-left)
n((X - 40e-6).^2 + (Y + 3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % Core 6 (top-right)
n((X + 40e-6).^2 + (Y - 3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;

% % Core 7 (top-left)
n((X - 40e-6).^2 + (Y - 3^0.5*40e-6).^2 < 2.5e-6^2) = n_core;
end

