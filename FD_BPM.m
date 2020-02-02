% Author: Madhu Veettikazhy
% Date:12 November 2018
% Gaussian beam (CW) electric field propagation from the
% through fibre using Finite Difference Beam Propagation Method
% in 2D
% Reference: An assessment of FD BPM method - Youngchul Chung
%
% Douglas-Gunn Alternating Direction Implicit approach is based on me690-lctr-nts.pdf Eqs. 3.21, 3.22, 3.23a and 3.23b.
% Implicit steps are performed using the Thomson algorithm (https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm)
% ***************************************************************************************************
format long
format compact

% Test of execution speed, 10 steps through a 2000x2000 grid, only the one
% mandatory update:
% 
% calctype 1: 657 seconds
% calctype 2: 48.7 seconds
% calctype 3: 1.56 seconds
% 
% That means method 3 is 421 times faster than method 1
% 
% GPU speedup comparison on 2020-02-02, with 1000 steps through a 2000x2000 grid, one update:
% calctype 3: 114 seconds
% calctype 4: 83.0 seconds when compiling with GCC, slightly slower when compiling with MSVC
% calctype 5: 12.2 seconds
% The speedup from calctype 3 to 5 is 9x
% Thus, the speedup from calctype 1 to 5 is ~3800x

%% General parameters
lambda = 1e-6;          % [m] Wavelength
w_0 = 3e-6;             % [m] Initial waist plane 1/e^2 radius of the gaussian beam
Lz = 1e-2;              % [m] z propagation distance

n_cladding = 1.45;
n_core = 1.46;
core_radius = 15e-6;

%% Resolution-related parameters
targetzstepsize = 1e-6; % [m] z step size to aim for
Lx_main = 50e-6;        % [m] x side length of main area
Ly_main = 50e-6;        % [m] y side length of main area
Nx_main = 200;          % x resolution of main area
Ny_main = 200;          % y resolution of main area

%% Solver-related parameters
calctype = 5; % 1: MATLAB vectorized method (full 2D Crank-Nicolson), 2: MATLAB sequential method (2x1D Douglas-Gunn-ADI Crank-Nicolson), 3: Same as 2 but in mex, 4: Same as 3 but with single-precision, 5: same as 4 but on GPU
useAllCPUs = true;

absorbertype = 4;         % 1: No absorber, 2: constant absorber, 3: linear absorber, 4: quadratic absorber
targetLx = 1.5*Lx_main;   % [m] Full area x side length, including absorber layer
targetLy = 1.5*Ly_main;   % [m] Full area y side length, including absorber layer
alpha = 10e1;             % [1/m] "Absorption coefficient", used for absorbertype 2
beta = 10e4;              % [1/m^2] "Absorption coefficient" per unit length distance out from edge of main area, used for absorbertype 3
gamma = 3e14;             % [1/m^3] "Absorption coefficient" per unit length distance out from edge of main area, squared, used for absorbertype 4

%% Visualization parameters
updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state. If set higher than nz (such as Inf), script will simply update every step.
colormax = 1e10;          % Maximum to use for the color scale in figure 4a

% expectedRayleighrange = pi*w_0^2/(lambda/n_core)
% expectedspotsize = w_0*sqrt(1+(Lz/expectedRayleighrange)^2)

%% Check for GPU compatibility if needed
if (calctype == 5)
  v = ver;
  if ~any(strcmp({v.Name},'Parallel Computing Toolbox'))
    error('You must have the Parallel Computing Toolbox installed to use GPU acceleration');
  end
  try
    GPUDev = gpuDevice;
    if(str2double(GPUDev.ComputeCapability) < 3)
      error('Your GPU is too old (CUDA compute capability < 3.0)')
    end
  catch
    error('No supported NVIDIA GPU found, or its driver is too old');
  end
end

%% Initialization of space and frequency grids
nz = round(Lz/targetzstepsize) % Number of z steps

dx = Lx_main/Nx_main;
dy = Ly_main/Ny_main;
dz = Lz/nz;

nx = round(targetLx/dx);
if rem(nx,2) ~= rem(Nx_main,2)
  nx = nx + 1; % Ensure that if Nx_main was set odd (to have a x slice at the center), nx will also be odd
end
ny = round(targetLy/dy);
if rem(ny,2) ~= rem(Ny_main,2)
  ny = ny + 1; % Ensure that if Ny_main was set odd (to have a y slice at the center), ny will also be odd
end
Lx = nx*dx;
Ly = ny*dy;

x = dx*(-(nx-1)/2:(nx-1)/2);
y = dy*(-(ny-1)/2:(ny-1)/2);
[X,Y] = ndgrid(x,y);

switch absorbertype
  case 1
    absorber = ones(nx,ny);
  case 2
    absorber = ones(nx,ny);
    absorber(abs(X) > Lx_main/2 | abs(Y) > Ly_main/2) = exp(-dz*alpha);
  case 3
    absorber = exp(-dz*max(0,max(abs(Y) - Ly_main/2,abs(X) - Lx_main/2))*beta);
  case 4
    absorber = exp(-dz*max(0,max(abs(Y) - Ly_main/2,abs(X) - Lx_main/2)).^2*gamma);
end
absorb_column = absorber(:);
figure(1);clf;
imagesc(x,y,absorber.');
axis xy
axis equal
xlim([-Lx/2 Lx/2]);
ylim([-Ly/2 Ly/2]);
xlabel('x [m]');
ylabel('y [m]');
colorbar;
title('Absorber illustration');

k_0 = 2*pi/lambda; % [m^-1] Wavenumber

%% Initialisation of optical fibre parameters
n_0 = n_core;
n_mat = ones(nx,ny)*n_cladding;
% n_mat(sqrt((X-Lx_main/4).^2+Y.^2) < core_radius) = n_core;
% n_mat(sqrt((X+Lx_main/4).^2+Y.^2) < core_radius) = n_core;
% n_mat(abs(X-Lx_main/4) < core_radius & abs(Y) < core_radius) = n_core;
% n_mat(abs(X+Lx_main/4) < core_radius & abs(Y) < core_radius) = n_core;
n_mat(sqrt(X.^2+Y.^2) < core_radius) = n_core;
% n_mat = max(n_cladding,n_core-(n_core-n_cladding)/core_radius^2*(X.^2+Y.^2));

figure(2);clf;
imagesc(x,y,n_mat.');
axis xy
axis equal
xlim([-Lx/2 Lx/2]);
ylim([-Ly/2 Ly/2]);
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title('Refractive index');

%% Beam initialization
% amplitude = exp(-((X-Lx_main/4).^2+Y.^2)/w_0^2) - exp(-((X+Lx_main/4).^2+Y.^2)/w_0^2); % Gaussian field amplitude
amplitude = exp(-((X-Lx_main/10).^2+Y.^2)/w_0^2); % Gaussian field amplitude
% amplitude = exp(-(X.^2+Y.^2)/w_0^2); % Gaussian field amplitude
phase = 8e5*Y;
E = amplitude.*exp(1i*phase); % Electric field
E = complex(E/sqrt(dx*dy*sum(abs(E(:)).^2)));
if(calctype == 4 || calctype == 5)
  E = complex(single(E));
end

%% Fibre propagation and plotting
updatezindices = round((1:min(nz,updates))/min(nz,updates)*nz);
z_updates = dz*updatezindices;
nextupdatenumber = 1;

switch calctype
  case 1
    delta_n_2 = n_mat.^2 - n_0^2; %delta_n^2 in the expression for FD BPM
    delta_n_2_column = delta_n_2(:); % Converting from NxN square matrix to N^2x1 column matrix
    E_column = E(:); % Converting 2D matrix into a column matrix
    absorber_column = absorber(:);
    N = nx*ny; % N*N - size of sparse matrices
    
    ax = dz/(2*dx^2);
    ay = dz/(2*dy^2);
    b =  dz*(1/dx^2+1/dy^2) - (k_0^2*dz/2).*delta_n_2_column + (2*1i*k_0*n_0);
    c = -dz*(1/dx^2+1/dy^2) + (k_0^2*dz/2).*delta_n_2_column + (2*1i*k_0*n_0);
    M_i_lhs = sparse(1:N,1:N,b(1:N),N,N);
    M_iMinus1 = sparse(2:N,1:N-1,ax,N,N);
    M_iPlus1 = sparse(1:N-1,2:N,ax,N,N);
    for i = nx:nx:N-nx % Dirichlet's boundary condition
      M_iMinus1(i+1,i) = 0;
      M_iPlus1(i,i+1) = 0;
    end
    
    M_jMinus1 = sparse(nx+1:N,1:N-nx,ay,N,N);
    M_jPlus1 = sparse(1:N-nx,nx+1:N,ay,N,N);
    M_lhs = -M_iMinus1 - M_iPlus1 + M_i_lhs - M_jMinus1 - M_jPlus1;
    
    M_i_rhs = sparse(1:N,1:N,c(1:N),N,N);
    M_rhs = M_iMinus1+M_iPlus1+M_i_rhs+M_jMinus1+M_jPlus1;
  case 2
    ax = dz/(4i*dx^2*k_0*n_0);
    ay = dz/(4i*dy^2*k_0*n_0);
    
    multiplier = exp(-1i*dz*k_0/2*(n_mat.^2 - n_0^2)/n_0).*absorber;
    d = zeros(nx,ny);
    b = zeros(max(nx,ny),1);
  case 3
    ax = dz/(4i*dx^2*k_0*n_0);
    ay = dz/(4i*dy^2*k_0*n_0);
    
    multiplier = complex(exp(-1i*dz*k_0/2*(n_mat.^2 - n_0^2)/n_0).*absorber);
    parameters = struct('multiplier',multiplier,'ax',ax,'ay',ay,'useAllCPUs',useAllCPUs);
  case {4,5}
    ax = single(dz/(4i*dx^2*k_0*n_0));
    ay = single(dz/(4i*dy^2*k_0*n_0));
    
    multiplier = complex(single(exp(-1i*dz*k_0/2*(n_mat.^2 - n_0^2)/n_0).*absorber));
    parameters = struct('multiplier',multiplier,'ax',ax,'ay',ay,'useAllCPUs',useAllCPUs);
end

%% Figure initialization
figure(3);clf;
h_line = semilogy(x,abs(E(:,round((ny-1)/2+1)).').^2);
h_title = title(['x intensity profile at z = 0 m (second moment radius = ' num2str(2*std(x,h_line.YData)) ' m)']);
line([-Lx_main/2 -Lx_main/2 NaN Lx_main/2 Lx_main/2],[ylim NaN ylim],'Linestyle','--','Color','k');
axis([min(x) max(x) 1e-9*colormax 10*colormax]);
xlabel('x [m]');
ylabel('Intensity [W/m^2]');

figure(4);clf;
subplot(2,1,1);
hold on;
h_im1 = imagesc(x,y,abs(E.').^2);
axis xy;
axis equal;
xlim([-Lx/2 Lx/2]);
ylim([-Ly/2 Ly/2]);
colorbar;
caxis('manual');
xlabel('x [m]');
ylabel('y [m]');
title('Intensity [W/m^2]');
if max(n_mat(:) > min(n_mat(:))); contour(X,Y,n_mat,(n_cladding+eps(n_cladding))*[1 1],'color','w','linestyle','--'); end
caxis([0 colormax]);
subplot(2,1,2);
hold on;
h_im2 = imagesc(x,y,angle(E.'));
axis xy;
axis equal;
xlim([-Lx/2 Lx/2]);
ylim([-Ly/2 Ly/2]);
colorbar;
caxis([-pi pi]);
if max(n_mat(:) > min(n_mat(:))); contour(X,Y,n_mat,(n_cladding+eps(n_cladding))*[1 1],'color','w','linestyle','--'); end
xlabel('x [m]');
ylabel('y [m]');
title('Phase [rad]');
colormap(gca,hsv/1.5);

powers = NaN(1,min(nz,updates)+1);
powers(1) = sum(dx*dy*abs(E(:)).^2);
figure(5);clf;
plot([0 z_updates],powers,'YDataSource','powers','linewidth',2);
xlim([0 Lz]);
xlabel('Propagation distance [m]');
ylabel('Relative power remaining');

tic
switch calctype
  case 1
    %% Vectorized method execution
    for zidx = 1:nz
      RHS = M_rhs*E_column;
      E_column = M_lhs\RHS;
      E_column = E_column.*absorber_column;
      if zidx == updatezindices(nextupdatenumber)
        E = reshape(E_column,[nx ny]);
        powers(nextupdatenumber+1) = dx*dy*sum(abs(E(:)).^2);
        h_line.YData = abs(E(:,round((ny-1)/2+1))).^2;
        h_title.String = ['x intensity profile at z = ' num2str(zidx*dz,'%.1e') ' m (second moment radius = ' num2str(2*std(x,h_line.YData),'%.6e') ' m)'];
        h_im1.CData = abs(E.').^2;
        h_im2.CData = angle(E.');
        refreshdata(5);
        drawnow;
        nextupdatenumber = nextupdatenumber + 1;
      end
    end
  case 2
    %% Sequential method execution
    for zidx = 1:nz
      %% Explicit part of step 1
      for iy = 1:ny
        for ix = 1:nx % construct d, the right hand side of the system of linear equations
          delta = 0;
          if ix ~= 1;  delta = delta +   ax*(E(ix-1,iy)-E(ix,iy)); end
          if ix ~= nx; delta = delta +   ax*(E(ix+1,iy)-E(ix,iy)); end
          if iy ~= 1;  delta = delta + 2*ay*(E(ix,iy-1)-E(ix,iy)); end
          if iy ~= ny; delta = delta + 2*ay*(E(ix,iy+1)-E(ix,iy)); end
          d(ix,iy) = E(ix,iy) + delta;
        end
      end
      
      %% Implicit part of step 1, Thomas algorithm along x, sweeps up from 2 to nx and then down from nx to 1
      for iy = 1:ny
        % Forward sweep, index 1 and 2:
        b(1) = 1 + ax;
        w = -ax/b(1);
        b(2) = 1 + 2*ax + w*ax;
        d(2,iy) = d(2,iy) - w*d(1,iy);
        % Forward sweep, indices 3 to nx-1:
        for ix = 3:nx-1
          w = -ax/b(ix-1);
          b(ix) = 1 + 2*ax;
          b(ix) = b(ix) + w*ax;
          d(ix,iy) = d(ix,iy) - w*d(ix-1,iy);
        end
        % Forward sweep, index nx:
        w = -ax/b(nx-1);
        b(nx) = 1 + ax;
        b(nx) = b(nx) + w*ax;
        d(nx,iy) = d(nx,iy) - w*d(nx-1,iy);
        
        % Back sweep, index nx:
        d(nx,iy) = d(nx,iy)/b(nx);
        % Back sweep, indices nx-1 to 2:
        for ix = nx-1:-1:2
          d(ix,iy) = (d(ix,iy) + ax*d(ix+1,iy))/b(ix);
        end
        % Back sweep, index 1:
        d(1,iy) = (d(1,iy) + ax*d(2,iy))/b(1);
      end
      
      %% Explicit part of step 2
      for iy = 1:ny
        for ix = 1:nx % construct d, the right hand side of the system of linear equations
          i = sub2ind([nx,ny],ix,iy);
          delta = 0;
          if iy ~= 1;  delta = delta - ay*(E(ix,iy-1)-E(ix,iy)); end
          if iy ~= ny; delta = delta - ay*(E(ix,iy+1)-E(ix,iy)); end
          d(ix,iy) = d(ix,iy) + delta;
        end
      end
      
      %% Implicit part of step 2, Thomas algorithm along y, sweeps up from 2 to ny and then down from ny to 1
      for ix = 1:nx
        % Forward sweep, index 2:
        b(1) = 1 + ay;
        w = -ay/b(1);
        b(2) = 1 + 2*ay + w*ay;
        d(ix,2) = d(ix,2) - w*d(ix,1);
        % Forward sweep, indices 3 to ny-1:
        for iy = 3:ny-1
          w = -ay/b(iy-1);
          b(iy) = 1 + 2*ay;
          b(iy) = b(iy) + w*ay;
          d(ix,iy) = d(ix,iy) - w*d(ix,iy-1);
        end
        % Forward sweep, index ny:
        w = -ay/b(ny-1);
        b(ny) = 1 + ay;
        b(ny) = b(ny) + w*ay;
        d(ix,ny) = d(ix,ny) - w*d(ix,ny-1);
        
        % Back sweep, index ny:
        E(ix,ny) = d(ix,ny)/b(ny);
        % Back sweep, indices ny-1 to 2:
        for iy = ny-1:-1:2
          E(ix,iy) = (d(ix,iy) + ay*E(ix,iy+1))/b(iy);
        end
        % Back sweep, index 1:
        E(ix,1) = (d(ix,1) + ay*E(ix,2))/b(1);
      end
      E = E.*multiplier;
      
      if zidx == updatezindices(nextupdatenumber)
        powers(nextupdatenumber+1) = dx*dy*sum(abs(E(:)).^2);
        h_line.YData = abs(E(:,round((ny-1)/2+1))).^2;
        h_title.String = ['x intensity profile at z = ' num2str(zidx*dz,'%.1e') ' m (second moment radius = ' num2str(2*std(x,h_line.YData),'%.6e') ' m)'];
        h_im1.CData = abs(E.').^2;
        h_im2.CData = angle(E.');
        refreshdata(5);
        drawnow;
        nextupdatenumber = nextupdatenumber + 1;
      end
    end
  case 3
    for updidx = 1:length(updatezindices)
      if updidx == 1
        parameters.nz = updatezindices(1);
      else
        parameters.nz = updatezindices(updidx) - updatezindices(updidx-1);
      end
      E = FDBPMpropagator(E,parameters);
      
      powers(updidx+1) = dx*dy*sum(abs(E(:)).^2);
      h_line.YData = abs(E(:,round((ny-1)/2+1))).^2;
      h_title.String = ['x intensity profile at z = ' num2str(updatezindices(updidx)*dz,'%.1e') ' m (second moment radius = ' num2str(2*std(x,h_line.YData),'%.6e') ' m)'];
      h_im1.CData = abs(E.').^2;
      h_im2.CData = angle(E.');
      refreshdata(5);
      drawnow;
    end
  case {4,5}
    for updidx = 1:length(updatezindices)
      if updidx == 1
        parameters.nz = updatezindices(1);
      else
        parameters.nz = updatezindices(updidx) - updatezindices(updidx-1);
      end
      if calctype == 4
        E = FDBPMpropagator_floats(E,parameters);
      else
        E = FDBPMpropagator_floats_CUDA(E,parameters);
      end
      
      powers(updidx+1) = dx*dy*sum(abs(E(:)).^2);
      h_line.YData = abs(E(:,round((ny-1)/2+1))).^2;
      h_title.String = ['x intensity profile at z = ' num2str(updatezindices(updidx)*dz,'%.1e') ' m (second moment radius = ' num2str(2*std(x,h_line.YData),'%.6e') ' m)'];
      h_im1.CData = abs(E.').^2;
      h_im2.CData = angle(E.');
      refreshdata(5);
      drawnow;
    end
%     figure(6);clf;imagesc(log10(abs(E.')));
end
toc

finalwidth = 2*std(x,h_line.YData)
finalpower = powers(end)