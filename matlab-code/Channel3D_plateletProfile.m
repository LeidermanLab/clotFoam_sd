% clear all
close all
format long

h_chan = 50e-3; % mm
w_chan = 150e-3; % mm

% In the inputParameters file
shear_y_low = 300; % 1/s low shear test (top/bottom of channel)
shear_z_low = 220; % Low shear test (front/back of channel)

shear_y_high = 1500; % 1/s High shear test (top/bottom of channel)
shear_z_high = 1100; % High shear test (front/back of chanel)


% % For the H channel in Link et al.
% shear_y = 1600; % 1/s High shear test (channel height direction)
% shear_z = 1300; % High shear test (channel width direction)
% h_chan = 60e-3; % mm
% w_chan = 100e-3; % mm

y0 = 0;
y_max = y0 + h_chan;
z0 = 0; 
z_max = z0 + w_chan;

positions = [
    0.095 0.25 0.35 0.6;
    0.5 0.25 0.35 0.6;
];

figure()
subplot(1,2,1)
inletProfile(shear_y_low,shear_z_low,h_chan,w_chan,y0,y_max,z0,z_max)

% clim([8.152616423454795e+04 7.844890633611609e+05])
xlabel('z (\mum)')
ylabel('y (\mum)')
set(gca, 'Position', positions(1, :));

subplot(1,2,2)
inletProfile(shear_y_high,shear_z_high,h_chan,w_chan,y0,y_max,z0,z_max)
colorbar()
xlabel('z (\mum)')
set(gca, 'Position', positions(2, :));
annotation(gcf,'textbox',...
    [0.94 0.86 0.05 0.1],...
    'String','P^{m,u} (plts/mm^3)',...
    'FontSize',16,...
    'FitBoxToText','on',...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'Rotation',270);

% Figure positioning
pos = [3 1 1100 230];
set(gcf, 'Position', pos);

function [] = inletProfile(shear_y,shear_z,h_chan,w_chan,y0,y_max,z0,z_max)

P0 = 2.5e5; % Conc. of platelets [plt/mm^3]

% Vessel radius in x and z directions
ry = h_chan / 2;
rz = w_chan / 2;

% --- uniform grid for OpenFOAM (uses cell centers)
Ny = 100; %18; % number of cells in y-direction
dy = h_chan / (Ny); % Cell height (magSf()/dz in OpenFoam)
y = y0 + dy/2: dy : y_max-dy/2; % face centers (Cf[2][i] in OpenFoam)

Nz = 300; %87; % number of cells in z-direction
dz = w_chan / (Nz); % Cell width (magSf()/dy in OpenFoam)
z = z0 + dz/2: dz : z_max-dz/2; % face centers (Cf[0][i] in OpenFoam)

[Z,Y] = meshgrid(z,y);

% Area of each cell on the boundary (magSf() in OpenFoam)
dA = dy*dz;

% To determine K, we observed that there is a linear relationship.
K_C = polyfit(500:500:1500,[202 330 455],2);
K_fxn = @(s) K_C(1)*s.^2 + K_C(2)*s + K_C(3);
Ky = K_fxn(shear_y);
Kz = K_fxn(shear_z);

% --- shape function ---
% c(y,z) = C0 [ 1 + K*Ry^{m-1}*(1-Ry)^{n-1} 
%               + K*Rz^{m-1}*(1-Rz)^{n-1}], 
% with m = 19, n = 2.
% Functions come from Eckstein 1991 (see Karins 2011 references)
S = zeros(Ny,Nz);
for i=1:Ny
for j=1:Nz
   Rz = abs(Z(i,j) - z0 - rz)/rz;
   Ry = abs(Y(i,j) - y0 - ry)/ry;
   S(i,j) = 1 + Ky * Ry^18 * (1-Ry) + Kz * Rz^18 * (1-Rz);
end
end

% Determine the max of s(y,z) at z = (z_max - z0)/2
Max_S = max(S(:,round(Nz/2)));
S = min(S,Max_S);
% for i=1:Ny
% for j=1:Nz
%    if (S(i,j) > Max_S)
%     S(i,j) = Max_S;
%    end
% end
% end

% Integrate the shape function over [y0, y_max] w/midpoint rule
intgrl_s = 0;
for i = 1:Ny
for j = 1:Nz
    intgrl_s = intgrl_s + dA*S(i,j);
end
end

% Scale for dimensional use
C0 = (w_chan*h_chan)/intgrl_s;
P=P0*C0*S;

% Display pertinent information
Ky = Ky
Kz = Kz
Max_S = Max_S
Max_P = max(max(P))
% peakToCenterRatio_y = max(P(round(Ny/2),:))/min(P(round(Ny/2),:)) % key as shear rate changes
% peakToCenterRatio_z = max(P(:,round(Nz/2)))/min(P(:,round(Nz/2),:)) % key as shear rate changes

pcolor(Z,Y,P)
shading interp
colormap('turbo')
% xlabel('z (\mum)')
% ylabel('y (\mum)')
% colorbar()
clim([8.152616423454795e+04 7.844890633611609e+05])
xlim([0 .15])
ylim([0 0.05])
xticks([0 0.05 0.1 0.15])
xticklabels({'0','50','100','150'})
yticks([0 0.05])
yticklabels({'0','50'})
ax = gca; 
ax.FontSize = 16; 
end

% % Calculate int_A P dA for constant P0
% intconst = P0 * (w_chan*h_chan)
% 
% % Check integral for new platelet profile (should = intconst)
% intprofile = 0;
% for k = 1:Ny
% for i = 1:Nz
%     intprofile = intprofile + dA*P(k,i);
% end
% end
% intprofile