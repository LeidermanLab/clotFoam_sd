clear all
close all
format long
hold on

% In the inputParameters file
shear_z = 1600; % 1/s
h_chan = 60e-3; % mm
w_chan = 100e-3; % mm
x0 = 0; %250e-3;
x_max = x0 + w_chan;
z0 = 0;
z_max = z0 + h_chan;
P0 = 2.5e5; % Conc. of platelets [plt/mm^3]

% Vessel radius in x and z directions
rx = w_chan / 2;
rz = h_chan / 2;

% --- uniform grid for OpenFOAM (uses cell centers)
Nx = 1000; %100; % number of cells in y-direction
dx = w_chan / (Nx); % Cell height (magSf()/dz in OpenFoam)
x = x0 + dx/2: dx : x_max-dx/2; % face centers (Cf[0][i] in OpenFoam)

Nz = 500; %50; % number of cells in y-direction
dz = h_chan / (Nz); % Cell height (magSf()/dx in OpenFoam)
z = z0 + dz/2: dz : z_max-dz/2; % face centers (Cf[2][i] in OpenFoam)

[X,Z] = meshgrid(x,z);

% Area of each cell on the boundary (magSf() in OpenFoam)
dA = dx*dz;

% To determine K, we observed that there is a linear relationship.  Using
% least squares, we define f_K(shear) = 0.253*shear + 76.0
A = [500 1; 1000 1; 1500 1]; % Matrix that multiplies slopeInt = [m b]'
Ks = [202 330 455]'; % from Karin's shape function
slopeInt = (A'*A) \ (A'*Ks);
f_K = @(x) slopeInt(1) * x + slopeInt(2);

K = f_K(shear_z);

% --- shape function 
% c(x,z) = C0 [ 1 + K*Rx^{m-1}*(1-Rx)^{n-1} 
%               + K*Rz^{m-1}*(1-Rz)^{n-1}], 
% with m = 19, n = 2.
% Functions come from Eckstein 1991 (see Karins 2011 references)
S = zeros(Nz,Nx);
for k=1:Nz
for i=1:Nx
   Rx = abs(X(k,i) - x0 - rx)/rx;
   Rz = abs(Z(k,i) - z0 - rz)/rz;
   S(k,i) = 1 + K * Rx^18 * (1-Rx) + K * Rz^18 * (1-Rz);
end
end

Max_S = max(S(:,round(Nx/2)));
for k=1:Nz
for i=1:Nx
   if (S(k,i) > Max_S)
    S(k,i) = Max_S;
   end
end
end

% Integrate the shape function over [y0, y_max] w/midpoint rule
intgrl_s = 0;
for k = 1:Nz
for i = 1:Nx
    intgrl_s = intgrl_s + dA*S(k,i);
end
end

% Scale for dimensional use
C0 = (w_chan*h_chan)/intgrl_s;
P=P0*C0*S;

% Display pertinent information

K = K
Max_S = Max_S
% peakToCenterRatio_x = max(P(round(Nz/2),:))/min(P(round(Nz/2),:)) % key as shear rate changes
peakToCenterRatio_z = max(P(:,round(Nx/2)))/min(P(:,round(Nx/2),:)) % key as shear rate changes

surf(X,Z,P)
shading interp
colormap('jet')

% Calculate int_A P dA for constant P0
intconst = P0 * (w_chan*h_chan)

% Check integral for new platelet profile (should = intconst)
intprofile = 0;
for k = 1:Nz
for i = 1:Nx
    intprofile = intprofile + dA*P(k,i);
end
end
intprofile