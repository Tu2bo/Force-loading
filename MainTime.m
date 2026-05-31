%% Main simulation script for the BR-integrated motor-clutch model
% This script implements the minimal mechanochemical model developed in:
%
% Xue R, Kang L, Chen Y, Yang H, Jiang H, Gong Z.
% "Force loading on molecular clutches governs the stability of cell lamellipodia."
% Proceedings of the National Academy of Sciences, 2026.
% DOI: https://doi.org/10.1073/pnas.2604349123
%
% This script performs:
% (i) Stochastic Monte Carlo simulations of the BR-integrated motor-clutch model;
% (ii) Mean-field ODE simulations for direct comparison;
% (iii) Quantitative analysis of the force-loading rate, cell spreading radius,
%       and force magnitude per bond.
%
% All model parameters and their physical meanings are provided in the main text
% and Supporting Information of the published article.
%
% Developed and tested in MATLAB R2024a.
% Last updated: May 2026

clc;
clear;

%% ---------------- Model parameters ----------------
kc    = 1;        % Clutch stiffness
kon   = 2.5;      % Association rate
koff0 = 0.5;      % Zero-force dissociation rate
Nc    = 500;      % Number of clutches
fb    = 2;        % Characteristic unbinding force

fm = [0.5 0.3 0.1] * Nc * fb;   % Total myosin force levels
r0 = 5e3;                       % Initial cell radius

% ---------------- Actin & membrane parameters ----------------
delta = 2.7;     % Intercalation gap
N     = 6000;    % Total actin filaments
M     = 100;     % Filaments per adhesion
kbT   = 4;       % Thermal energy
sigma = 0.02;    % Membrane tension
vu    = 150;     % Unloaded retrograde flow speed

kmem = 4*pi*sigma*M/N;          % Effective membrane stiffness
rp   = N*kbT/(4*pi*sigma*delta);% Characteristic polymerization length

% ---------------- Substrate & kinetic parameters ----------------
ks     = 7e7;                  % Substrate stiffness
vbeta  = fb*koff0/kc;          % Characteristic velocity scale
alpha  = 1.35;                 % Force sensitivity parameter

% ---------------- Simulation settings ----------------
Tfinal = 1e4;
y01    = [0 0 r0];              % Initial conditions for mean-field model
Tspan  = [0 Tfinal];
options = odeset('RelTol',1e-6,'AbsTol',1e-8);


%%
parfor i = 1:length(fm)
    % Monte Carlo simulation
    res(i,1) = MC_BR_motor_clutch_model (kc, ks, vu, kon, koff0, fb, Nc, fm(i), kmem, rp, r0, Tfinal);

    % Mean-field simulation
     [T{i},y{i}] = ode15s(@(t,y) MF_BR_motor_clutch_model(t, y, kc, ks, vu, kon, koff0, fb, Nc, fm(i), kmem, rp, vbeta, alpha), Tspan, y01,options);
end

%% ---------------- Loading rate ----------------
INTV = 2e4;

for i = 1:length(res)
    figure; hold on;

    LR_mc = ks*kc*res(i).vf ./ (ks + res(i).Nbt*kc);
    LR_mc(res(i).Nbt < 1) = inf;

    scatter(res(i,1).T(1:INTV:end),LR_mc(1:INTV:end),'filled','SizeData',12,'MarkerFaceColor','#257FC9')

    vf_mf = vu*(1 - (y{i}(:,1) - kmem*y{i}(:,3))/fm(i));
    LR_mf = ks*kc*vf_mf ./ (ks + y{i}(:,2)*kc);
    LR_mf(y{i}(:,2) < 1) = inf;

    plot(T{i},LR_mf,'LineWidth',1.5,'color',"#78B64E")

    Lcrit = lambertw(vu./vbeta) * fb * kon / alpha;
    yline(Lcrit, 'k');

    set(gca,'FontSize',8,'Color','none','TickLength',[0.02 0.02],'YMinorTick','on','XMinorTick','on','box','on','Layer','top','Yscale','log')
    set(gcf,"Units","centimeters","Position",[30,10,5.5,2.5])
    xlim([0 Tfinal])
    ylim([1e-1 1e3])
    yticks([1e-1 1e1 1e3])
end

%% ---------------- Cell radius ----------------
for i = 1:length(res)
    figure; hold on;

    scatter(res(i,1).T(1:INTV:end),res(i,1).r(1:INTV:end)*1e-3,'filled','SizeData',12,'MarkerFaceColor','#257FC9')
    plot(T{i},y{i}(:,3)*1e-3,'LineWidth',1.5,'color',"#78B64E")

    set(gca,'FontSize',8,'Color','none',...
        'TickLength',[0.02 0.02],'YMinorTick','on','XMinorTick','on','box','on','Layer','top')
    set(gcf,"Units","centimeters","Position",[30,10,5.5,2.5])
    xlim([0 Tfinal])
    ylim([0 120])
    yticks(0:40:120)
end

%% ---------------- Force magnitude  ----------------
for i = 1:length(res)
    figure;hold on;

    scatter(res(i,1).T(1:INTV:end),res(i,1).fadh(1:INTV:end)./res(i,1).Nbt(1:INTV:end),'filled','SizeData',12,'MarkerFaceColor','#257FC9')
    plot(T{i},y{i}(:,1)./y{i}(:,2),'LineWidth',1.5,'color',"#78B64E")

    Fcritical = (1+lambertw(kon/koff0/exp(1)))/alpha*fb;
    line([0 Tfinal],[Fcritical Fcritical],"color","k")

    set(gca,'FontSize',8,'Color','none',...
        'TickLength',[0.02 0.02],'YMinorTick','on','XMinorTick','on','box','on','Layer','top')
    set(gcf,"Units","centimeters","Position",[30,10,5.5,2.5])
    xlim([0 Tfinal])
    ylim([0 6])
    yticks(0:2:6)
end

