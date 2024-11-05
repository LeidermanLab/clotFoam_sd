% Phenomenological Functions (linear interpolation)

% Adhesion parameters
kAdhVwfOn_low = 2.9891e-08; % on-rate of adhesion via vWF [mm^3 / s]
kAdhVwfOn_high = 8.7181e-08; % on-rate of adhesion via vWF [mm^3 / s]
kAdhVwfOff_low = 0.4698; % off-rate of adhesion via vWF [1 / s]
kAdhVwfOff_high = 1.6; % off-rate of adhesion via vWF [1 / s]

% Cohesion parameters
kCohVwfOnPmax_low = 12000; % on-rate of cohesion via vWF [1/s]
kCohVwfOnPmax_high = 20000; % on-rate of cohesion via vWF [1/s]
kCohVwfOff_low = 374.813; % off-rate of cohesion via vWF [1/s]
kCohVwfOff_high = 500; % off-rate of cohesion via vWF [1/s]

% Shear rates
lowShear = 300;
highShear = 1500;
maxShear = 2500;

f = @(S, S_L, K_L, S_H, K_H) K_L * (S-S_H)./(S_L-S_H) + ...
                                  K_H*(S-S_L)./(S_H-S_L);  
f_min = @(S, S_L, K_L, S_H, K_H, S_M) min(K_L * (S-S_H)./(S_L-S_H) + ...
                                  K_H*(S-S_L)./(S_H-S_L),...
                                  K_L * (S_M-S_H)./(S_L-S_H) + ...
                                  K_H*(S_M-S_L)./(S_H-S_L));  


shear = linspace(0,5000,1000);

figure()
makePlot(f_min,shear,lowShear,kAdhVwfOn_low,highShear,kAdhVwfOn_high,maxShear,'k_{{adh}^+}^{vWF}(\gamma)  [mm^3 / s]')
function [] = makePlot(f,S,S_L,K_L,S_H,K_H,S_M,dispName)
K = f(S, S_L, K_L, S_H, K_H,S_M);
plot(S,K,'-k',linewidth=2)
hold on
plot([S_L,S_H],[K_L,K_H],'ko',linewidth=2,MarkerFaceColor='k')
xlabel('Shear Rate [1/s]')
ylabel(dispName)
set(gca,FontSize=16)
end
