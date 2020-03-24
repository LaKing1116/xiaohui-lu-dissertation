% Authors: Madhu Veetikazhy and Anders K. Hansen
% DTU Health and DTU Fotonik
% 
% Electric field propagation through fibre using Finite Difference Beam
% Propagation Method in 2D Reference: An assessment of FD BPM method -
% Youngchul Chung
%
% Douglas-Gunn Alternating Direction Implicit approach is based on
% me690-lctr-nts.pdf Eqs. 3.21, 3.22, 3.23a and 3.23b. Implicit steps are
% performed using the Thomson algorithm
% (https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm)
% ***********************************************************************

format long
format compact

FileName = 'singleMode_w0';  % File name for the saved video and data files
videoName = [FileName '.avi']; 
 
saveVideo = false;  % To save the field intensity and phase profiles at different transverse planes
saveData = false;  % To save the required variables from the simulation result

%% USER DEFINED General parameters
clear Lz taperScaling twistRate shapeTypes shapeParameters shapeRIs

lambda = 980e-9;          % [m] Wavelength
w_0 = 2.35e-6;             % [m] Initial waist plane 1/e^2 radius of the gaussian beam

n_cladding = 1.45;
n_core = 1.46;

% Lz, taperScaling, twistRate, shapeTypes, shapeParameters and shapeRIs are
% cell arrays in which each element corresponds to a fiber segment in the
% simulation. If an entry of shapeTypes is an empty array, it means the
% shapes should simply carry over from the previous segment. Otherwise, new
% shapes are defined, emulating a fiber splice.

Lz{1} = 0.15e-3; % [m] z propagation distances, one for each segment
taperScaling{1} = 1; % Specifies how much the refractive index profile of the last z slice should be scaled relative to the first z slice, linearly scaling in between
twistRate{1} = 0; % Specifies how rapidly the fiber twists, measured in radians per metre
shapeTypes{1} = [1]; % Shape types for each segment. An empty array in a cell means that the previous shapes carry over. Shape types are 1: Circular step-index disk, 2: Antialiased circular step-index disk, 3: Parabolic graded index disk
shapeParameters{1} = [0; % x values
                      0; % y values
                      2e-6]; % r values
shapeRIs{1} = [n_core]; % Refractive indices to use for the shapes

% Lz{2} = 5e-3;
% taperScaling{2} = 0.15;
% twistRate{2} = 2*pi/Lz{2};
% shapeTypes{2} = [];
% shapeParameters{2} = [];
% shapeRIs{2} = [];
% 
% Lz{3} = 2e-3;
% taperScaling{3} = 0.15;
% twistRate{3} = 0;
% shapeTypes{3} = [];
% shapeParameters{3} = [];
% shapeRIs{3} = [];
% 
% Lz{4} = 3e-3;
% taperScaling{4} = 0.15;
% twistRate{4} = 0;
% shapeTypes{4} = [3];
% shapeParameters{4} = [0.15*(2e-6)
%                       0.15*(12e-6)
%                       0.15*(10e-6)];
% shapeRIs{4} = [1.465];

Eparameters = {w_0};    % Cell array of parameters that the E field initialization function (defined at the end of this file) will need

%% USER DEFINED Resolution-related parameters
targetzstepsize = 0.5e-6; % [m] z step size to aim for
Lx_main = 25e-6;        % [m] x side length of main area
Ly_main = 25e-6;        % [m] y side length of main area
Nx_main = 400;          % x resolution of main area
Ny_main = 400;          % y resolution of main area

%% USER DEFINED Solver-related parameters
useAllCPUs = true;
useGPU = false;

n_0 = n_core;

targetLx = 1.5*Lx_main;   % [m] Full area x side length, including absorber layer
targetLy = 1.5*Ly_main;   % [m] Full area y side length, including absorber layer
alpha = 3e14;             % [1/m^3] "Absorption coefficient" per unit length distance out from edge of main area, squared

%% USER DEFINED Visualization parameters
updatesTotal = 300;            % Number of times to update plot. Must be at least 1, showing the final state.
colormax = 1;          % Maximum to use for the color scale in figure 3a
downsampleImages = false; % Due to a weird MATLAB bug, MATLAB may crash when having created imagesc (or image) plots with dimensions larger than roughly 2500x2500 and then calling mex functions repeatedly. This flag will enable downsampling to 500x500 of all data before plotting, hopefully avoiding the issue.
displayScaling = 6;  % Zooms in on figures 1 & 3a,b. Set to 2 for no zooming.  
if saveVideo
    video = VideoWriter(videoName);  
    open(video);
end






%% YOU SHOULDN'T NORMALLY NEED TO MODIFY ANYTHING BELOW THIS LINE, EXCEPT THE E INITALIZATION FUNCTION AT THE END OF THE FILE
%% Check for GPU compatibility if needed
if useGPU
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

%% Calculate updates for each segment
LzTotal = sum(cell2mat(Lz));
updates = cell(1,numel(Lz));
updates{1} = max(1,round(Lz{1}/LzTotal*updatesTotal));
LzRemaining = LzTotal - Lz{1};
nUpdatesRemaining = max(1,updatesTotal - updates{1});
for iSeg=2:numel(Lz) % Segment index
  updates{iSeg} = max(1,round(Lz{iSeg}/LzRemaining*nUpdatesRemaining));
  LzRemaining = LzRemaining - Lz{iSeg};
  nUpdatesRemaining = nUpdatesRemaining - updates{iSeg};
end
updatesCumSum = cumsum(cell2mat(updates));

%% Initialization of space and frequency grids
dx = Lx_main/Nx_main;
dy = Ly_main/Ny_main;

Nx = round(targetLx/dx);
if rem(Nx,2) ~= rem(Nx_main,2)
  Nx = Nx + 1; % Ensure that if Nx_main was set odd (to have a x slice at the center), Nx will also be odd
end
Ny = round(targetLy/dy);
if rem(Ny,2) ~= rem(Ny_main,2)
  Ny = Ny + 1; % Ensure that if Ny_main was set odd (to have a y slice at the center), Ny will also be odd
end
Lx = Nx*dx;
Ly = Ny*dy;

x = dx*(-(Nx-1)/2:(Nx-1)/2);
y = dy*(-(Ny-1)/2:(Ny-1)/2);
[X,Y] = ndgrid(x,y);

if downsampleImages
  if Nx>500
    ix_plot = round(linspace(1,Nx,500));
  else
    ix_plot = 1:Nx;
  end
  x_plot = x(ix_plot);

  if Ny>500
    iy_plot = round(linspace(1,Ny,500));
  else
    iy_plot = 1:Ny;
  end
  y_plot = y(iy_plot);
end

k_0 = 2*pi/lambda; % [m^-1] Wavenumber

%% Beam initialization
E = calcInitialE(X,Y,Eparameters); % Call function to initialize E field
E = complex(single(E/sqrt(max(abs(E(:)).^2)))); % Normalize and force to be complex single precision
% E = complex(single(E/sqrt(dx*dy*sum(abs(E(:)).^2)))); % Normalize and force to be complex single precision
E_0 = E;  % For initial intensity, phase, and power values

%% Figure initialization
figure(1);clf;
figure1_Settings = [];
screenSize = get(0,'MonitorPositions');
figure1_Settings.monitorNumber = 1; % Main monitor number
figure1_Settings.size.figureHeight = screenSize(figure1_Settings.monitorNumber,4); % Figure height (pixels)
figure1_Settings.size.figureWidth = screenSize(figure1_Settings.monitorNumber,3); % Figure width (pixels)
set(gcf, 'Position',  [0, 0, figure1_Settings.size.figureWidth, figure1_Settings.size.figureHeight]);

subplot(2,2,1)
if downsampleImages
  h_im1 = imagesc(x_plot,y_plot,zeros(min(500,Ny),min(500,Nx),'single'));
else
  h_im1 = imagesc(x,y,zeros(Ny,Nx,'single'));
end
axis xy
axis equal
xlim([-Lx/displayScaling Lx/displayScaling]);
ylim([-Ly/displayScaling Ly/displayScaling]);
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title('Refractive index');

powers = NaN(1,updatesCumSum(end)+1);
powers(1) = sum(abs(E(:)).^2)/sum(abs(E_0(:)).^2);
% powers(1) = sum(dx*dy*abs(E(:)).^2);
zUpdates = zeros(1,updatesCumSum(end)+1);
subplot(2,2,2);
plot(zUpdates,powers,'XDataSource','zUpdates','YDataSource','powers','linewidth',2);
xlim([0 LzTotal]);
xlabel('Propagation distance [m]');
ylabel('Relative power remaining');

subplot(2,2,3);
hold on;
box on;
if downsampleImages
  h_im3a = imagesc(x_plot,y_plot,abs(E(ix_plot,iy_plot).').^2);
else
  h_im3a = imagesc(x,y,abs(E.').^2);
end
axis xy;
axis equal;
xlim([-Lx/displayScaling Lx/displayScaling]);
ylim([-Ly/displayScaling Ly/displayScaling]);
colorbar;
caxis('manual');
xlabel('x [m]');
ylabel('y [m]');
title('Intensity [W/m^2]');
% if max(n_mat(:) > min(n_mat(:))); contour(X,Y,n_mat,(n_cladding+eps(n_cladding))*[1 1],'color','w','linestyle','--'); end
line([-Lx_main Lx_main Lx_main -Lx_main -Lx_main]/2,[Ly_main Ly_main -Ly_main -Ly_main Ly_main]/2,'color','r','linestyle','--');
caxis([0 colormax]);
colormap(gca,GPBGYRcolormap);

h_axis3b = subplot(2,2,4);
hold on;
box on;
if downsampleImages
  h_im3b = imagesc(x_plot,y_plot,angle(E(ix_plot,iy_plot).'/E(Nx/2+1+round(shapeParameters{1}(1)/dx),Ny/2+1+round(shapeParameters{1}(2)/dy)))); 
else
  h_im3b = imagesc(x,y,angle(E.'/E(Nx/2+1+round(shapeParameters{1}(1)/dx),Ny/2+1+round(shapeParameters{1}(2)/dy)))); 
end
h_axis3b.Color = 0.7*[1 1 1];  % To set the color corresponding to phase outside the cores where there is no field at all
axis xy;
axis equal;
xlim([-Lx/displayScaling Lx/displayScaling]);
ylim([-Ly/displayScaling Ly/displayScaling]);
colorbar;
caxis([-pi pi]);
h_im3b.AlphaData = max(0,(1+log10(abs(E.'/max(E(:))).^2)/3));  %Logarithmic transparency in displaying phase outside cores
% if max(n_mat(:) > min(n_mat(:))); contour(X,Y,n_mat,(n_cladding+eps(n_cladding))*[1 1],'color','w','linestyle','--'); end
line([-Lx_main Lx_main Lx_main -Lx_main -Lx_main]/2,[Ly_main Ly_main -Ly_main -Ly_main Ly_main]/2,'color','r','linestyle','--');
xlabel('x [m]');
ylabel('y [m]');
title('Phase [rad]');
colormap(gca,hsv/1.5);

drawnow;

if saveVideo
    frame = getframe(gcf);  %Get the frames 
    writeVideo(video,frame);  %Stitch the frames to form a video and save
end

%% Loop over the segments
NzTotal = 0;
tic;
for iSeg = 1:numel(Lz) % Segment index
  %% Calculate the segment-by-segment taper scaling factor to pass into the mex function
  if iSeg == 1
    segTaperScaling = taperScaling{iSeg};
  else
    segTaperScaling = taperScaling{iSeg}/taperScaling{iSeg-1};
  end
  %% Either redefine the shapes defining the RI profile (if this element of shapeTypes is non-empty) or rescale and rotate the previous one
  if ~isempty(shapeTypes{iSeg})
    segShapeTypes = shapeTypes{iSeg};
    segShapeParameters = shapeParameters{iSeg};
    segShapeRIs = shapeRIs{iSeg};
  else
    if iSeg == 2
      scaleFactor = taperScaling{iSeg-1};
    else
      scaleFactor = taperScaling{iSeg-1}/taperScaling{iSeg-2};
    end
    oldSegShapeParameters = segShapeParameters;
    segShapeParameters(:,1) = scaleFactor*(cos(twistRate{iSeg-1}*Lz{iSeg-1})*oldSegShapeParameters(:,1) - sin(twistRate{iSeg-1}*Lz{iSeg-1})*oldSegShapeParameters(:,2));
    segShapeParameters(:,2) = scaleFactor*(sin(twistRate{iSeg-1}*Lz{iSeg-1})*oldSegShapeParameters(:,1) + cos(twistRate{iSeg-1}*Lz{iSeg-1})*oldSegShapeParameters(:,2));
    segShapeParameters(:,3) = scaleFactor*oldSegShapeParameters(:,3);
  end
  
  %% Calculate z step size and positions
  Nz = max(updates{iSeg},round(Lz{iSeg}/targetzstepsize)); % Number of z steps in this segment
  NzTotal = NzTotal + Nz;
  dz = Lz{iSeg}/Nz;

  zUpdateIdxs = round((1:updates{iSeg})/updates{iSeg}*Nz); % Update indices relative to start of this segment
  if iSeg == 1
    zUpdates(2:updatesCumSum(iSeg)+1) = dz*zUpdateIdxs;
  else
    zUpdates(updatesCumSum(iSeg-1)+2:updatesCumSum(iSeg)+1) = dz*zUpdateIdxs + zUpdates(updatesCumSum(iSeg-1)+1);
  end

  %% Calculate proportionality factors for use in the mex function
  ax = dz/(4i*dx^2*k_0*n_0);
  ay = dz/(4i*dy^2*k_0*n_0);
  d = -dz*k_0/(2*n_0); % E = E*multiplier*exp(1i*d*(n^2-n_0^2))

  %% Define the multiplier
  absorber = exp(-dz*max(0,max(abs(Y) - Ly_main/2,abs(X) - Lx_main/2)).^2*alpha);
  multiplier = absorber; % This could also include a phase gradient due to bending

  %% Load variables into a parameters struct and start looping, one iteration per update
  parameters = struct('dx',single(dx),'dy',single(dy),'taperPerStep',single((1-segTaperScaling)/Nz),'twistPerStep',single(twistRate{iSeg}*Lz{iSeg}/Nz),...
    'shapeTypes',uint8(segShapeTypes),'shapeParameters',single(segShapeParameters),'shapeRIs',single(segShapeRIs),'n_cladding',single(n_cladding),'multiplier',complex(single(multiplier)),...
    'd',single(d),'n_0',single(n_0),'ax',single(ax),'ay',single(ay),'useAllCPUs',useAllCPUs);

  parameters.iz_start = int32(0); % z index of the first z position to step from for the first call to FDBPMpropagator, in C indexing (starting from 0)
  parameters.iz_end = int32(zUpdateIdxs(1)); % last z index to step into for the first call to FDBPMpropagator, in C indexing (starting from 0)
  for updidx = 1:length(zUpdateIdxs)
    if updidx > 1
      parameters.iz_start = int32(zUpdateIdxs(updidx-1));
      parameters.iz_end   = int32(zUpdateIdxs(updidx));
    end
    if useGPU
      [E,n] = FDBPMpropagator_CUDA(E,parameters);
    else
      [E,n] = FDBPMpropagator(E,parameters);
    end

    %% Update figure contents
    if downsampleImages
      h_im1.CData = n(ix_plot,iy_plot).'; % Refractive index at this update
      h_im3a.CData = abs(E(ix_plot,iy_plot).').^2; % Intensity at this update
      h_im3b.CData = angle(E(ix_plot,iy_plot).'/E(Nx/2+1+round(shapeParameters{1}(1)/dx),Ny/2+1+round(shapeParameters{1}(2)/dy))); % Phase at this update
    else
      h_im1.CData = n.'; % Refractive index at this update
      h_im3a.CData = abs(E.').^2; % Intensity at this update
      h_im3b.CData = angle(E.'/E(Nx/2+1+round(shapeParameters{1}(1)/dx),Ny/2+1+round(shapeParameters{1}(2)/dy))); % Phase at this update
    end
    
    if iSeg == 1
      powers(updidx+1) = sum(abs(E(:)).^2)/sum(abs(E_0(:)).^2);
%       powers(updidx+1) = dx*dy*sum(abs(E(:)).^2);
    else
      powers(updatesCumSum(iSeg-1) + updidx + 1) = sum(abs(E(:)).^2)/sum(abs(E_0(:)).^2);
%       powers(updatesCumSum(iSeg-1) + updidx + 1) = dx*dy*sum(abs(E(:)).^2);
    end
    refreshdata(1);
    drawnow;
    if saveVideo
        frame = getframe(1); 
        writeVideo(video,frame); 
    end
  end
end
toc

if saveVideo
	close(video);
end

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
% amplitude = exp(-((X-Lx_main/4).^2+Y.^2)/w_0^2) - exp(-((X+Lx_main/4).^2+Y.^2)/w_0^2); % Gaussian field amplitude
% amplitude = exp(-((X-Lx_main/10).^2+Y.^2)/w_0^2); % Gaussian field amplitude
% amplitude = exp(-(X.^2+Y.^2)/w_0^2); % Gaussian field amplitude
w_0 = Eparameters{1};
amplitude1 = exp(-(X.^2+Y.^2)/w_0^2);
% amplitude2 = 2*exp(-((X+12e-6).^2+(Y+7e-6).^2)/w_0^2);
phase1 = zeros(size(X));
% phase2 = 8e5*Y;
E = amplitude1.*exp(1i*phase1);% + amplitude2.*exp(1i*phase2); % Electric field
end