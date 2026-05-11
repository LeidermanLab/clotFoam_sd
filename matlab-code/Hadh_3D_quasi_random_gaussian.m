clear
close all

% Print c++ code for xn and zn?
matToCpp = 0;

% Diameter of a platelet [mm]
Pdiam = 3e-3;          

% Spatial Parameters for H domain
x0 = 0;             % Lower bound of x dimension [mm]
x_max = 350e-3;     % Upper bound of x dimension [mm]
y0 = 0;             % Lower bound of y dimension [mm]
y_max = 60e-3;      % Upper bound of y dimension [mm]
z0 = 0;             % Lower bound of z dimension [mm]
z_max = 80e-3;      % Upper bound of z dimension [mm]

% Injury dimensions
l_inj = 150e-3;     % Length of the injury in x-direction [mm]
w_inj = 60e-3;     % Width of the injury in the z-direction [mm]
distFromEdge = Pdiam/4; % The closest a Gaussian can get to edge of injury

% Number of grid cells in each direction
Nx = 100;            % number of cells in x-direction
Ny = 50;            % number of cells in y-direction
Nz = 100;            % number of cells in x-direction

% % Spatial Parameters for 3D rectangular channel
% x0 = 0;             % Lower bound of x dimension [mm]
% x_max = 160e-3;     % Upper bound of x dimension [mm]
% y0 = 0;             % Lower bound of y dimension [mm]
% y_max = 50e-3;      % Upper bound of y dimension [mm]
% z0 = 0;             % Lower bound of z dimension [mm]
% z_max = 150e-3;     % Upper bound of z dimension [mm]
% 
% % Injury dimensions
% l_inj = 100e-3;     % Length of the injury in x-direction [mm]
% w_inj = 100e-3;     % Width of the injury in the z-direction [mm]
% distFromEdge = Pdiam; % The closest a Gaussian can get to edge of injury
% 
% % Number of grid cells in each direction
% Nx = 68;            % number of cells in x-direction
% Ny = 50;            % number of cells in y-direction
% Nz = 87;            % number of cells in x-direction

% Gaussian Parameters 
sxz = 0.6;
sigma_x = sxz*Pdiam;  % Standard deviation in x (+/- 1 mu) [mm]
sigma_y = 3*Pdiam;       % Standard deviation in x (+/- 1 mu) [mm]
sigma_z = sxz*Pdiam;  % Standard deviation in z (+/- 1 mu) [mm]


% Halton Point Parameters
Np = 175;          
base_x = 2;         % Base for the Halton sequence in x dimension (max 5)
base_z = 3;         % Base for the Halton sequence in z dimension (max 6)

% nth Gaussian function gn(z,t) centered at (xn,zn)
gn = @(x,y,z,xn,zn) exp( -(x - xn).^2/(2*sigma_x^2) ... 
                         - y.^2/(2*sigma_y^2) ...
                         - (z-zn).^2/(2*sigma_z^2));

% Generate the Halton points in 2D
xn = zeros(Np,1);
zn = zeros(Np,1);
for i = 1:Np
        xn(i) = halton_points(i, base_x);
        zn(i) = halton_points(i, base_z);
end

% Scale Halton points for where injury is located
scale = @(a,b,x) a + (b-a)*x;
l_start = (x_max - l_inj)/2;
l_end = (x_max + l_inj)/2;
w_start = (z_max - w_inj)/2;
w_end = (z_max + w_inj)/2;
xn = scale(l_start+distFromEdge,l_end-distFromEdge,xn);
zn = scale(w_start+distFromEdge,w_end-distFromEdge,zn);

% Uniform grid for OpenFOAM (uses cell centers)
dx = (x_max - x0) / Nx; 
x = x0 + dx/2: dx : x_max-dx/2; % face centers (Cf[0][i] in OpenFoam)
dy = (y_max - y0) / Ny; 
y = y0 + dy/2: dy : y_max-dy/2; % face centers (Cf[1][i] in OpenFoam)
dz = (z_max - z0) / Nz; 
z = z0 + dz/2: dz : z_max-dz/2; % face centers (Cf[2][i] in OpenFoam)

% Calculate G in 3D
[X,Y,Z] = ndgrid(x,y,z);
G = 0*X;
for n = 1:Np
    % G = G + gn(X,Y,Z,xn(n),zn(n));
    G = min(1,G + gn(X,Y,Z,xn(n),zn(n)));
end
for j = 1:Ny
    if y(j) >= Pdiam
        G(:,j,:) = 0;
    end
end
maxG = max(max(max(G)));
G = G/maxG;

% Cpp code:
% forAll(CC,cellI) 
% if x_start < x < x_end, y < Pdiam, z_start < z < z_end
% Calculate G at particular x,y,z, then set Hadh

% Write necessary data to file
if matToCpp == 1
    fname = 'haltonPoints.dat';
    fid = fopen(fname,'wt');
    matNumToC(Np,'Np',fid);
    matNumToC(maxG,'maxG',fid);
    matNumToC(sigma_x^2,'sigma_x2',fid);
    matNumToC(sigma_y^2,'sigma_y2',fid);
    matNumToC(sigma_z^2,'sigma_z2',fid);
    matVecToC(xn,'xn',fid);
    matVecToC(zn,'zn',fid);
    fclose(fid);
end

% Calculate the average Hadh value
format long
avg_Hadh = avgHadh(G,dx,dy,dz,l_start,l_end,0,Pdiam,w_start,w_end)

filePath = '/Users/dmontgomery/Dropbox/Documents/College/Mines/Research/WorkWithDave/Weekly Notes/20240311-haltonPointsHadh/figs';

% Plot Halton Points
figure()
scatter(xn, zn, 'ko');
hold on
rectangle('Position',[l_start w_start l_inj w_inj])
xlabel('x [mm]');
ylabel('z [mm]');
title([num2str(Np),' Halton Points on x-z Plane']);
set(gca,'FontSize',16)
xlim([0 x_max])
ylim([0 z_max])
grid on;
% saveas(gcf,fullfile(filePath,['Hadh_',num2str(Np),'_HaltonPts']),'fig')
% saveas(gcf,fullfile(filePath,['Hadh_',num2str(Np),'_HaltonPts']),'png')

% Plot G at y = y(y_indx)
y_indx = 3;
figure()
Gxz = reshape(G(:,y_indx,:),Nx,Nz)';
[Xxz,Zxz] = meshgrid(x,z);
pcolor(Xxz,Zxz,Gxz)
hold on 
rectangle('Position',[l_start w_start l_inj w_inj],LineWidth=3)
shading flat
xlabel('x [mm]');
ylabel('z [mm]');
colorbar()
title(['Hadh(x,y,z) at y = ',num2str(1e3*y(y_indx)),' \mum']);
set(gca,'FontSize',16)
shading interp
% colormap('jet');
% saveas(gcf,fullfile(filePath,['Hadh_',num2str(Np),'_xz']),'fig')
% saveas(gcf,fullfile(filePath,['Hadh_',num2str(Np),'_xz']),'png')

% % Plot G at z = z(z_indx)
% z_indx = 43;
% figure()
% Gxy = reshape(G(:,:,z_indx),Nx,Ny)';
% [Xxy,Yxy] = meshgrid(x,y);
% pcolor(Xxy,Yxy,Gxy)
% hold on
% rectangle('Position',[l_start 0 l_inj Pdiam],LineWidth=3)
% shading flat
% xlabel('x [mm]');
% ylabel('y [mm]');
% colorbar()
% title(['Hadh(x,y,z), with z = ',num2str(1e3*z(z_indx)),' \mum']);
% set(gca,'FontSize',16)

% -----------------------------------------------------------------------
% Functions used throughout calculations above
% -----------------------------------------------------------------------
function points = halton_points(n, base)
    % Parameters:
    % n: Sequence index
    % base: Base for the Halton sequence
    
    points = 0;
    f = 1;
    
    while n > 0
        f = f / base;
        points = points + f * mod(n, base);
        n = floor(n / base);
    end
end

function [] = matVecToC(var,var_name,fid)
    [m,~] = size(var);
    if m > 1
        var = var';
    end
    precision=16;
    N=length(var);
    fmt=['double %s[%d]={' repmat('%d,',1,numel(var)-1) '%d};\n'];
    % c_code=sprintf(fmt,var_name,N,var)
    
    % fname = ['plt_Pmu_',var_name,'.dat'];
    % fid = fopen(fname,'wt');
    fprintf(fid,fmt,var_name,N,var);
    % fclose(fid);
end

function [] = matNumToC(var,var_name,fid)
    if mod(var,1) == 0
        fmt='int %s = %d;\n';
    else
        fmt='scalar %s = %d;\n';
    end
    fprintf(fid,fmt,var_name,var);
    % fname = ['plt_Pmu_',var_name,'.dat'];
    % fid = fopen(fname,'wt');
    % fprintf(fid,fmt,var_name,N,var);
    % fclose(fid);
end

function H_bar = avgHadh(G,dx,dy,dz,x_start,x_end,y_start,y_end,z_start,z_end)
    
    index = @(x,d) max(round(x/d + 0.5),1);

    i_start = index(x_start,dx); i_end = index(x_end,dx);
    j_start = index(y_start,dy); j_end = index(y_end,dy);
    k_start = index(z_start,dz); k_end = index(z_end,dz);

    H = G(i_start:i_end,j_start:j_end,k_start);
    for k = k_start+1:k_end
        H = [H; G(i_start:i_end, j_start:j_end, k)];
    end
    
    H_bar = mean(mean(H));

end
