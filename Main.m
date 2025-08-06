% Main function for Monte Carlo simulation and main field theory
% The BR-integrated motor-clutch model
% Dimensionless 
% Ruihao Xue 2025
clc
clear
Nc = 500;                                                          % Number of clutches
R0 = 111.8;                                                       % Initial cell radius
V = 150;                                                            % Maximal retrograde flow speed
Vp0 = 6.71;                                                       % Maximal ploymerization speed
Fp = 790.85;                                                     % Characteristic ploymerization force
Ron = 5;                                                           % On-rate
Kmem = 1.87e-4;                                             % Membrane stiffness
alpha = 1.35;                                                    % Correction coefficient
Fm = [0.15 0.3 0.35 0.5];                                  % Myosin force
K = 1;                                                               % Substrate stiffness
Tfinal = 2e3;
y01 = [0 0 R0];
Tspan1 = [0 Tfinal];
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
%%
parfor i = 1:length(Fm)
    res(i,1) = SingleSpringEntDim(Ron, K, V, Vp0, Fm(i), Kmem, Fp, Nc, R0, Tfinal); %Stochastic simulation 
    [T{i},y{i}] = ode15s(@(t,y) MFSingleSpringEntDim(t, y, K, V, Vp0, Fp, Ron, Fm(i), Kmem, alpha), Tspan1, y01,options); %Mean-field theory 
end

%%
for i = 1 : length(Fm)
    Hs = fsolve(@(H) Fm(i)*(1-H)-Kmem*Fp*log(H)-Ron*lambertw(V*H)./(Ron+exp(lambertw(V*H)))/alpha,1);
    Fs(i) =  lambertw(V*Hs)/alpha;

    figure (i); %F/P force-loading magnitude
    hold on
    plot(res(i,1).T,res(i,1).F,'LineWidth',1.5) % Stochastic simulation
    plot(T{i},y{i}(:,1),'LineWidth',1.5)% Mean-field theory

    set(gca,'FontSize',10,'Color','none',...
        'TickLength',[0.02 0.02],'Layer','top','YMinorTick','on','XMinorTick','on','box','on')
    set(gcf,"Units","centimeters","Position",[30,10,4,3])
    xlim([0 Tfinal])
    ylim([0, 0.8])
    yticks(0:0.4:0.8)

    figure (i+length(Fm)); %KVH force-loading rate
    hold on
    Hm =  1-(res(i,1).F-Kmem*res(i,1).R)/Fm(i);
    plot(res(i,1).T,K*V*Hm,'LineWidth',1.5)
    H = 1-(y{i}(:,1)-Kmem*y{i}(:,3))/Fm(i);
    plot(T{i},K*V*H,'LineWidth',1.5)
    set(gca,'FontSize',10,'Color','none','Yscale','log',...
        'TickLength',[0.02 0.02],'Layer','top','YMinorTick','on','XMinorTick','on','box','on')
    set(gcf,"Units","centimeters","Position",[30,10,4,3])
    ylim([1 3e2])
    yticks([1e-1 1e0 1e1 1e2])
end

%% Cell radius
figure (17)
hold on
for i = length(Fm):-1:1
    plot(res(i,1).T,res(i,1).R,'LineWidth',1.5)  % Stochastic simulation 
    % plot(T{i},y{i}(:,3),'LineWidth',1.5) % Mean-field theory
end
set(gca,'FontSize',10,'Color','none',...
    'TickLength',[0.02 0.02],'Layer','top','YMinorTick','on','XMinorTick','on','box','on')
set(gcf,"Units","centimeters","Position",[30,10,6,3])
xlim([0 Tfinal])
ylim([0 2500])
yticks(0:1200:2500)

