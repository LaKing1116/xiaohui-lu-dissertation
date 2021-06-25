function P = findModes(P,nModes,singleCoreModes,sortByLoss,plotModes)
if isfield(P,'n_cladding')
  error('Error: n_cladding has been renamed n_background');
end
if isfield(P,'n') && isfield(P,'shapes')
  error('Error: You must specify exactly one of the fields "shapes" and "n"');
end
if isfield(P,'n')
  if ~isfield(P.n,'n') && ~isfield(P.n,'func')
    error('If you specify "n", it must contain either the field "func" or "n"');
  end
  if isfield(P.n,'n') && ~(isfield(P.n,'Lx') && isfield(P.n,'Ly'))
    error('You must specify the side lengths Lx and Ly if you provide an array (2D or 3D) for the refractive index');
  end
  if isfield(P.n,'func') && nargin(P.n.func) == 5 && ~(isfield(P.n,'Lx') && isfield(P.n,'Ly') && isfield(P.n,'Nx') && isfield(P.n,'Ny') && isfield(P.n,'Nz'))
    error('You must specify the side lengths Lx and Ly as well as the refractive index array resolutions Nx, Ny and Nz if you provide a 3D refractive index function');
  end
end
if isfield(P,'shapes') && size(P.shapes,2) == 5
  if any(P.shapes(:,4) == 4) || any(P.shapes(:,4) == 5)
    error('Since you have a GRIN lens, you must define the gradient constant g in the shapes array');
  else
    P.shapes(:,6) = NaN;
  end
end
if ~isfield(P,'nParameters')
  P.nParameters = {};
end
if ~isfield(P,'rho_e')
  P.rho_e = 0.22;
end
if ~isfield(P,'bendingRoC')
  P.bendingRoC = Inf;
end
if ~isfield(P,'bendDirection')
  P.bendDirection = 0;
end
if ~isfield(P,'Intensity_colormap')
  P.Intensity_colormap = 1;
end
if ~isfield(P,'Phase_colormap')
  P.Phase_colormap = 2;
end

k_0 = 2*pi/P.lambda;  % [m^-1] Wavenumber
dx = P.Lx_main/P.Nx_main;
dy = P.Ly_main/P.Ny_main;

targetLx = P.padfactor*P.Lx_main;
targetLy = P.padfactor*P.Ly_main;

Nx = round(targetLx/dx);
if rem(Nx,2) ~= rem(P.Nx_main,2)
  Nx = Nx + 1; % Ensure that if Nx_main was set odd (to have a x slice at the center), Nx will also be odd
end
Ny = round(targetLy/dy);
if rem(Ny,2) ~= rem(P.Ny_main,2)
  Ny = Ny + 1; % Ensure that if Ny_main was set odd (to have a y slice at the center), Ny will also be odd
end
Lx = Nx*dx;
Ly = Ny*dy;

x = dx*(-(Nx-1)/2:(Nx-1)/2);
y = dy*(-(Ny-1)/2:(Ny-1)/2);
[X,Y] = ndgrid(x,y);

V = [];
D = [];

if ~isfield(P,'shapes')
  shapesToInclude_2Darray = 1;
elseif singleCoreModes
  shapesToInclude_2Darray = (1:size(P.shapes,1)).';
else
  shapesToInclude_2Darray = 1:size(P.shapes,1);
end

fprintf('Finding modes...\n');
tic
for iModeFinderRun = 1:size(shapesToInclude_2Darray,1)
  n = P.n_background*ones(Nx,Ny);
  if isfield(P,'shapes')
    for iShape = shapesToInclude_2Darray(iModeFinderRun,:)
      switch P.shapes(iShape,4)
        case 1
          n((X-P.shapes(iShape,1)).^2 + (Y-P.shapes(iShape,2)).^2 < P.shapes(iShape,3)^2) = P.shapes(iShape,5);
        case 2
          delta = max(dx,dy);
          r_diff = sqrt((X-P.shapes(iShape,1)).^2 + (Y-P.shapes(iShape,2)).^2) - P.shapes(iShape,3) + delta/2;
          n(r_diff < delta) = min(r_diff(r_diff < delta)/delta*(P.n_background - P.shapes(iShape,5)) + P.shapes(iShape,5),P.shapes(iShape,5));
        case 3
          r_ratio_sqr = ((X-P.shapes(iShape,1)).^2 + (Y-P.shapes(iShape,2)).^2)/P.shapes(iShape,3)^2;
          n(r_ratio_sqr < 1) = r_ratio_sqr(r_ratio_sqr < 1)*(P.n_background - P.shapes(iShape,5)) + P.shapes(iShape,5);
        case 4
          r_ratio_sqr = ((X-P.shapes(iShape,1)).^2 + (Y-P.shapes(iShape,2)).^2)/P.shapes(iShape,3)^2;
          r_abs = sqrt((X-P.shapes(iShape,1)).^2 + (Y-P.shapes(iShape,2)).^2);
          n(r_ratio_sqr < 1) = 2*P.shapes(iShape,5)*exp(P.shapes(iShape,6)*r_abs(r_ratio_sqr < 1))./(exp(2*P.shapes(iShape,6)*r_abs(r_ratio_sqr < 1)) + 1);
        case 5
          r_ratio_sqr = (Y-P.shapes(iShape,2)).^2/P.shapes(iShape,3)^2;
          r_abs = Y - P.shapes(iShape,2);
          n(r_ratio_sqr < 1) = 2*P.shapes(iShape,5)*exp(P.shapes(iShape,6)*r_abs(r_ratio_sqr < 1))./(exp(2*P.shapes(iShape,6)*r_abs(r_ratio_sqr < 1)) + 1);
      end
    end
  elseif isfield(P.n,'func') % Otherwise if P.n has a function field
    if nargin(P.n.func) == 4 % It's 2D
      n = single(P.n.func(X,Y,P.n_background,P.nParameters));
    else % It's 3D. We evaluate it on the user requested grid but with only z = 0 and then interpolate out to the full grid
      d_n(1) = P.n.Lx/P.n.Nx;
      d_n(2) = P.n.Ly/P.n.Ny;
      d_n(3) = P.Lz/(P.n.Nz-1);
      x_n = d_n(1)*(-(P.n.Nx-1)/2:(P.n.Nx-1)/2);
      y_n = d_n(2)*(-(P.n.Ny-1)/2:(P.n.Ny-1)/2);
      z_n = 0;
      [X_n,Y_n,Z_n] = ndgrid(single(x_n),single(y_n),single(z_n));
      n = single(P.n.func(X_n,Y_n,Z_n,P.n_background,P.nParameters));
      clear X_n Y_n Z_n
      n = interpn(x_n,y_n,n,x,y.','linear',P.n_background);
    end
  else % Otherwise a P.n.n array must have been specified
    [Nx_source,Ny_source,Nz_source] = size(P.n.n);
    if Nz_source == 1 % If P.n.n is 2D, interpolate it to the simulation grid
      dx_source = P.n.Lx/Nx_source;
      dy_source = P.n.Ly/Ny_source;
      x_source = dx_source*(-(Nx_source-1)/2:(Nx_source-1)/2);
      y_source = dy_source*(-(Ny_source-1)/2:(Ny_source-1)/2);
      n = interpn(x_source,y_source,P.n.n,x,y.','linear',P.n_background);
    else % Otherwise it's 3D so we take the first slice and interpolate to the siulation grid
      d_n(1) = P.n.Lx/Nx_source;
      d_n(2) = P.n.Ly/Ny_source;
      x_n = d_n(1)*(-(P.n.Nx-1)/2:(P.n.Nx-1)/2);
      y_n = d_n(2)*(-(P.n.Ny-1)/2:(P.n.Ny-1)/2);
      n = interpn(x_n,y_n,P.n.n(:,:,1),x,y.','linear',P.n_background);
    end
  end

  n = double(n); % Needed for sparse operations
  
  anycomplex = ~isreal(n);
  
  n_bend = real(n).*(1-(real(n).^2.*(X*cosd(P.bendDirection) + Y*sind(P.bendDirection))*P.rho_e/(2*P.bendingRoC))).*exp((X*cosd(P.bendDirection) + Y*sind(P.bendDirection))/P.bendingRoC);

  delta_n_2 = real(n_bend).^2 - P.n_0^2;                              %delta_n^2 in the expression for FD BPM

  dz = 1e-10;
  absorber = exp(-dz*(max(0,max(abs(Y) - P.Ly_main/2,abs(X) - P.Lx_main/2)).^2*P.alpha + 2*pi*imag(n)/P.lambda)); % First part is edge absorber, second is absorption from the imaginary part of the refractive indicres
  ax = 1.00001*dz/(dx^2*2i*k_0*P.n_0);
  % ax = dz/(dx^2*2i*k_0*P.n_0);
  ay = dz/(dy^2*2i*k_0*P.n_0);


  N = Nx*Ny;                                                            %N*N - size of sparse matrices
  M_rhs = sparse(1:N,1:N,absorber(1:N) + delta_n_2(1:N)*dz*k_0/(2i*P.n_0),N,N) + ...
    sparse(1:N-1,2:N,[repmat([repmat(ax,1,Nx-1) 0],1,Ny-1) repmat(ax,1,Nx-1)],N,N) + ...
    sparse(2:N,1:N-1,[repmat([repmat(ax,1,Nx-1) 0],1,Ny-1) repmat(ax,1,Nx-1)],N,N) + ...
    sparse(1:N-Nx,1+Nx:N,ay,N,N) + ...
    sparse(1+Nx:N,1:N-Nx,ay,N,N);
  M_rhs(1:N+1:N*N) = M_rhs(1:N+1:N*N) - repmat([ax repmat(2*ax,1,Nx-2) ax],1,Ny);
  M_rhs(1:N+1:N*N) = M_rhs(1:N+1:N*N) - [repmat(ay,1,Nx) repmat(2*ay,1,Nx*(Ny-2)) repmat(ay,1,Nx)];
  absorberdiag = sparse(1:N,1:N,absorber(1:N),N,N);
  M_rhs = M_rhs*absorberdiag;

  [Vrun,Drun] = eigs(M_rhs,ceil(nModes/size(shapesToInclude_2Darray,1)),1,'Display',false,'SubspaceDimension',nModes*10);
  V = [V , Vrun];
  D = [D , diag(Drun).'];
end
fprintf('\b Done, %.1f seconds elapsed\n',toc);

if sortByLoss
  [~,sortedidxs] = sort(real(D),'descend');
else
  [~,sortedidxs] = sort(imag(D),'ascend');
end

% for iMode = nModes:-1:1
for iMode = 1:nModes
  P.modes(iMode).Lx = Lx;
  P.modes(iMode).Ly = Ly;
  E = reshape(V(:,sortedidxs(iMode)),[Nx Ny]);
  P.modes(iMode).field = E.*exp(-1i*angle(max(E(:)))); % Shift phase arbitrarily so that the resulting modes are (nearly) real
  kappa = (1-real(D(sortedidxs(iMode))))/dz*P.lambda/(4*pi);
  realn = sqrt(P.n_0^2 - 2*P.n_0*imag(D(sortedidxs(iMode))/(dz*k_0)));
  if anycomplex
    P.modes(iMode).neff = realn + kappa*1i;
  else
    P.modes(iMode).neff = realn; % If the user did not specify any complex refractive indices, then the only contribution to kappa would be from the edge absorber, and showing a complex neff would just confuse people and not be very physically meaningful
  end
  if isfield(P,'shapes') && (isinf(P.bendingRoC) && (size(P.shapes,1) == 1 || singleCoreModes))
    [~,iMax] = max(E(:));
    xMax = X(iMax);
    yMax = Y(iMax);
    theta = atan2(yMax - P.shapes(1,2),xMax - P.shapes(1,1));
    radialE = interp2(X.',Y.',E.',P.shapes(1,1) + linspace(0,max(Lx,Ly),1000)*cos(theta),P.shapes(1,2) + linspace(0,max(Lx,Ly),1000)*sin(theta));
    radialEpruned = radialE(abs(radialE) > 0.01*max(abs(radialE)));
    m = sum(abs(diff(angle(radialEpruned) > 0))) + 1;
    
    R = sqrt((xMax - P.shapes(1,1))^2 + (yMax - P.shapes(1,2))^2);
    azimuthalE = interp2(X.',Y.',E.',P.shapes(1,1) + R*cos(theta + linspace(0,2*pi,1000)),P.shapes(1,2) + R*sin(theta + linspace(0,2*pi,1000)));
    azimuthalEpruned = azimuthalE(abs(azimuthalE) > 0.01*max(abs(azimuthalE)));
    l = sum(abs(diff(angle(azimuthalEpruned) > 0)))/2;
    
    if l > 0
      Emaxmirrored = interp2(X.',Y.',E.',xMax,2*P.shapes(1,2) - yMax);
      if real(E(iMax)/Emaxmirrored) < 0
        parity = 'o';
      else
        parity = 'e';
      end
    else
      parity = '';
    end
    P.modes(iMode).label = ['Mode ' num2str(iMode) ', LP' num2str(l) num2str(m) parity];
  else
    P.modes(iMode).label = ['Mode ' num2str(iMode)];
  end
  if plotModes
    h_f = figure(100+iMode);
    h_f.WindowStyle = 'docked';
    subplot(1,2,1);
    imagesc(x,y,abs(E.').^2);
    if isfield(P,'shapes')
      for iShape = 1:size(P.shapes,1)
        if P.shapes(iShape,4) <=3
          line(P.shapes(iShape,1) + P.shapes(iShape,3)*cos(linspace(0,2*pi,100)),P.shapes(iShape,2) + P.shapes(iShape,3)*sin(linspace(0,2*pi,100)),'color','w','linestyle','--');
        end
      end
    end
    axis equal; axis tight; axis xy;
    setColormap(gca,P.Intensity_colormap);
    subplot(1,2,2);
    maxE0 = max(E(:));
    imagesc(x,y,angle(E.'/maxE0),'AlphaData',max(0,(1+log10(abs(E.'/maxE0).^2)/3)));
    if isfield(P,'shapes')
      for iShape = 1:size(P.shapes,1)
        if P.shapes(iShape,4) <=3
          line(P.shapes(iShape,1) + P.shapes(iShape,3)*cos(linspace(0,2*pi,100)),P.shapes(iShape,2) + P.shapes(iShape,3)*sin(linspace(0,2*pi,100)),'color','w','linestyle','--');
        end
      end
    end
    set(gca,'Color',0.7*[1 1 1]);  % To set the color corresponding to phase outside the cores where there is no field at all
    caxis([-pi pi]);
    axis equal; axis tight; axis xy;
    setColormap(gca,P.Phase_colormap);
    if anycomplex
      sgtitle({[P.modes(iMode).label ', n_{eff} = ' num2str(realn,'%.6g') ' + ' num2str(kappa,'%.3g') 'i'],['rough loss estimate: ' num2str(-log(real(D(sortedidxs(iMode))))/dz,'%.3g') ' m^{-1} (' num2str((-10*log10(exp(-1)))*(-log(real(D(sortedidxs(iMode))))/dz),'%.3g') ' dB/m)']});
    else
      sgtitle({[P.modes(iMode).label ', n_{eff} = ' num2str(realn,'%.6g')],['rough loss estimate: ' num2str(-log(real(D(sortedidxs(iMode))))/dz,'%.3g') ' m^{-1} (' num2str((-10*log10(exp(-1)))*(-log(real(D(sortedidxs(iMode))))/dz),'%.3g') ' dB/m)']});
    end
  end
end
drawnow;
end

function setColormap(gca,colormapType)
    switch colormapType
     case 1
        colormap(gca,GPBGYRcolormap);
     case 2
        colormap(gca,hsv/1.5);
     case 3
        colormap(gca,parula);
     case 4
        colormap(gca,gray);
     case 5
        colormap(gca,cividisColormap);
    end
end