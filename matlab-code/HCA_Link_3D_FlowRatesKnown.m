% This script solves the Hydraulic Circuit from Danes and Leiderman for use
% in simulating the results in the Link et al. paper which had many
% mistakes in the HCA.

clear all

% Density of blood [g/cm^3]
rho = 1;

% Prescribed flow rates [uL/min]
Q1 = 5.5;
Q3 = 10;

% Prescribed outlet pressure [Pa]
Pcw0 = 0;

% Pressure drop across the injury [Pa]
Pdrop = 2.80851e2;

% Lengths of the bleeding chip [micrometers]
W_legs = 100; % Width of channels (x-dir)
L_legs = 85; % (y-dir)
H_legs = 60; % (z-dir)
W_inj = 20; % (y-dir)
L_inj = 165; % Effective length of injury channel (x-dir)
H_inj = 60; % Height of the injury channel (z-dir)
L_up = 4915; % Length of the channel above comp. domain (y-dir)
% Viscosities by region [Poise]
mu_b = 3.6e-2; 
mu_up_wash = 1e-2; 
mu_inj = 2.67e-2;
mu_down_wash = 2.3e-2;

% --- Begin Calculations Don't Edit! -------------------------------------
% Convert lengths from mu m to m
W_legs = convert_from_mum_to_m(W_legs);
L_legs = convert_from_mum_to_m(L_legs);
H_legs = convert_from_mum_to_m(H_legs);
L_up = convert_from_mum_to_m(L_up);
W_inj = convert_from_mum_to_m(W_inj);
L_inj = convert_from_mum_to_m(L_inj);
H_inj = convert_from_mum_to_m(H_inj);

% Convert viscosity from Poise to Pa*s
mu_b = convert_from_Poise_to_Pa_s(mu_b);
mu_up_wash = convert_from_Poise_to_Pa_s(mu_up_wash);
mu_inj = convert_from_Poise_to_Pa_s(mu_inj);
mu_down_wash = convert_from_Poise_to_Pa_s(mu_down_wash);

% Calculate all of the resistances [Pa*s / m^3]
a = @(L,w,h) 12*(1 - 192*h/(pi^5*w)*tanh(pi*w/(2*h)))^(-1);
calc_R = @(L,w,h,mu) a(L,w,h)*mu*L/(w*h^3);

R1 = calc_R(L_legs,W_legs,H_legs,mu_b);
R2 = calc_R(L_legs,W_legs,H_legs,mu_b);
R3 = calc_R(L_legs,W_legs,H_legs,mu_up_wash);
R4 = calc_R(L_legs,W_legs,H_legs,mu_down_wash);
Rm = calc_R(L_inj,W_inj,H_inj,mu_inj);
Rub = calc_R(L_up,W_legs,H_legs,mu_b);
Ruw = calc_R(L_up,W_legs,H_legs,mu_up_wash);

% Convert flowrates from muL/min to m^3/s
Q1 = convert_from_uLPerMin_to_m3PerS(Q1);
Q3 = convert_from_uLPerMin_to_m3PerS(Q3);

% Calculate Qm [mm^3/s]
Qm = Pdrop / Rm;

% Build the matrix A
%      Q2    Q4    Pb   Pcb  Pcb0    P1    Pw   Pcw    P2      
A = [   1,    0,    0,    0,    0,    0,    0,    0,    0;
        0,    1,    0,    0,    0,    0,    0,    0,    0;
        0,    0,    1,   -1,    0,    0,    0,    0,    0;
        0,    0,    0,    1,    0,   -1,    0,    0,    0;
      -R2,    0,    0,    0,   -1,    1,    0,    0,    0;
        0,    0,    0,    0,    0,    0,    1,   -1,    0;
        0,    0,    0,    0,    0,    0,    0,    1,   -1;
        0,  -R4,    0,    0,    0,    0,    0,    0,    1;
        0,    0,    0,    0,    0,    1,    0,    0,   -1];


% RHS vector b
b = [Q1-Qm; Q3+Qm; Q1*Rub; Q1*R1; 0; Q3*Ruw; Q3*R3; Pcw0; Qm*Rm];

% Solve for unknowns
x = A\b;

% Various vectors for display
flowRateNames = ['Q1';'Q2';'Q3';'Q4';'QM'];
flowRates = [Q1; x(1); Q3; x(2); Qm];
resistanceNames = ['Rub';'Ruw';'R1 '; 'R2 '; 'R3 '; 'R4 '; 'Rm '];
R = [Rub; Ruw; R1; R2; R3; R4; Rm];
pressureNames = ['Pb   ';'Pw   ';'Pcb  ';'Pcb0 ';'Pcw  ';'Pcw0 ';'P1   ';'P2   ';'P1-P2'];
P = [     x(3);   x(7);   x(4);   x(5);   x(8);   Pcw0;    x(6); x(9); x(6) - x(9)];
% Pdrop_calc = P(7) - P(8)
% Pdrop

disp(' ')
disp('------------------------------------------')
disp('--- Hydraulic Circuit Analysis Results ---')
disp('------------------------------------------')

% format shortE
disp('Calculated Restances')
disp(table(resistanceNames,R,...
     'VariableNames',{'Label','Resistance (Pa s / m^3)'}))

disp('Flow Rates')
disp(table(flowRateNames,flowRates*6e10,flowRates,flowRates*1e9,...
     'VariableNames',{'Label','Flow Rate (uL / min)','Flow Rate (m^3 / s )', 'Flow Rate (mm^3 / s )'}))

rho = convert_from_gPerCm3_to_kgPerM3(rho);
kinPressure = P/rho; 

disp('Pressure ')
disp(table(pressureNames,P,kinPressure,convert_from_m2_to_mm2(kinPressure),...
     'VariableNames',{'Label','Pressure (Pa)', 'Kin. Pressure (m^2 / s^2)','Kin. Pressure (mm^2 / s^2)'}))
disp(['Pdrop = ', num2str(Pdrop)])

% disp('Kinematic Pressure ')



function lnew = convert_from_mum_to_m(L)
    % From mu m to m 
    lnew = 1e-6*L;
end

function lnew = convert_from_m_to_mm(L)
    % From m to mm
    lnew = L * 1e-2;
end

function lnew2 = convert_from_m2_to_mm2(L2)
    lnew2 = convert_from_m_to_mm(convert_from_m_to_mm(L2));
    lnew2 = 1e6*L2;
end


function mu_new = convert_from_Poise_to_Pa_s(mu)
    % Convert from Poise to Pa*s
    mu_new = 0.1*mu;
end


function Qnew = convert_from_uLPerMin_to_m3PerS(Q)
    % From mu m to m 
    Qnew = 1/6*1e-10*Q;
end

function rhoNew = convert_from_gPerCm3_to_kgPerM3(rho)
    rhoNew = rho*1e3;
end