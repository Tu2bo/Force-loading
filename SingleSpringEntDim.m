%Ruihao Xue; 2025
%The BR-integrated motor-clutch model
%Dimensionless
function results = SingleSpringEntDim(Ron , K, V, Vp0, Fm, Kmem, Fp, Nc, R0, Tfinal)
% definition of basic parameters
Pi = false(Nc,1);
Ron = Ron*ones(Nc,1);
Fci = zeros(Nc,1);
Uf = zeros(Nc,1);
Time = zeros(1000,1);
F = zeros(1000,1);
Pbt= zeros(1000,1);
Vf= zeros(1000,1);
Vp= zeros(1000,1);
R= R0*ones(1000,1);
% solving governing equations based on Monte Carlo method
n = 2;
T = 0;

while T<Tfinal
    % Determine the clutch event
    Roffi = exp(Fci);                                  % slip bond
    Ra = not(Pi).*Ron+Pi.*Roffi;
    Rtot = sum(Ra);
    Rnum = rand(2,1);
    dT = -log(Rnum(1))/Rtot;                                                        % Time step size
    dT = max(dT,1e-8);
    cum_Ra = cumsum(Ra/Rtot);                                                   % Cumulative Probability Distribution
    ind = sum(Rnum(2) > cum_Ra) + 1;                                        % Vectorized Search
    Pi(ind) = 1-Pi(ind);                                                                    % Changing the state
    Pb = sum(Pi)/Nc;
    H = 1- (F(n-1)-Kmem*R(n-1))/Fm;
    % Calculate the deformation of clutches and the subtrate
    Ufb = Pi.*(V*H*dT+Uf);
    Us = sum(Ufb)/(Nc*(K+Pb));

    % Calculate the clutch force and retrograde flow
    Uf = not(Pi).*Us +Ufb;
    Fci = Uf-Us;
    T = T+dT;
    %Uf is clutch displacement: for unbound clutches Uf=Us; for bounded Uf=Ufb
    F(n) = sum(Fci)/Nc;
    Pbt(n) = Pb;
    Time(n) = T;
    Vf(n) = Vp0*H;
    Vp(n) = Vp0*exp(-R(n-1)/Fp);
    R(n) = Vp0*(exp(-R(n-1)/Fp)-H)*dT+R(n-1);
    % R(n) = Vp0*(1-H)*dT+R(n-1);  %For the conventional motor-clutch model

    n = n+1;
end
results.T = Time;
results.F = F;
results.Pb = Pbt;
results.R = R;
results.Vf = Vf;
results.Vp = Vp;
end
