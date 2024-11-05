% close all
clear

% Hill functions for ADP
ADP_shift = 1e-5;
% ADP_ast = [5e-5 7.5e-5 1e-4 1.5e-4 2e-3];
ADP_ast = [1e-4 1.5e-4 2e-3];
maxADP = 0.001; % [nanomol/mm^3]
n = 1; % degree of Hill function


% convert to nM by multiplying by 10^6
ADP_ast = ADP_ast*1e6; % [nM] 2000 muM in Leiderman 2011 
maxADP = maxADP*1e6; % [nM]
s = ADP_shift*1e6;

ADP = linspace(0,maxADP,1000);
A = @(c,c_ast) max(0,(c-s).^n) ./ (c_ast.^n + max(0,(c-s).^n));

lineW = 2;
figure()
hold on
set(gca,'ColorOrderIndex',1)
for i = 1:length(ADP_ast)-1
plot(ADP,0.34*A(ADP,ADP_ast(i)),'-','DisplayName',['[ADP]^* = ',num2str(ADP_ast(i)),' nM'],LineWidth=lineW)
end

% s = 0;
% A = @(c,c_ast) max(0,(c-s).^n) ./ (c_ast.^n + max(0,(c-s).^n));

plot(ADP,20*0.34*A(ADP,ADP_ast(end)),'-','DisplayName',['[ADP]^* = ',num2str(ADP_ast(end)),' nM'],LineWidth=lineW)
xlabel('Concentration [nM]')
legend()
% title(['n = ', num2str(n)])
