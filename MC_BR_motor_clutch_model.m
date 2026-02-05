%
%Stochastic Monte Carlo simulations for the BR-integrated motor–clutch model;
%
function results = MC_BR_motor_clutch_model(kc, ks, vu, kon, koff0, fb, Nc, Fm, kmem, rp, R0, Tfinal)
% definition of basic parameters
Nb = 0;
Pbi = false(Nc,1);
koni = kon*ones(Nc,1);
fc = zeros(Nc,1);
uc = zeros(Nc,1);
vf = vu*ones(1000,1);                                                   % Retrograde flow speed (unit nm/s)
vp = vu*ones(1000,1);                                                   % (unit nm/s)
r = R0*ones(100,1);                                                % Cell spreading radius (unit nm)
T = zeros(1000,1);
fadh = zeros(1000,1);
fmeant = zeros(1000,1);
Nbt= zeros(1000,1);
% solving governing equations based on Monte Carlo method
n = 2;
t = 0;
while t<Tfinal
    fadh(n-1) = sum(fc);          % Total traction force (unit pN)
    Nbt(n-1) = Nb;
    fmeant(n-1) = fadh(n-1)/(Nbt(n-1)+eps);
    T(n-1) = t;

    % Determine the clutch event
    koffi = koff0*exp(fc./fb);                                  % slip bond
    ra = not(Pbi).*koni+Pbi.*koffi;
    rtot = sum(ra);
    rnum = rand(2,1);
    dt = -log(rnum(1))/rtot;                                                        % Time step size
    cum_ra = cumsum(ra/rtot);                                                   % Cumulative Probability Distribution
    ind = sum(rnum(2) > cum_ra) + 1;                                        % Vectorized Search
    Pbi(ind) = 1-Pbi(ind);                                                                % Changing the state
    Nb = sum(Pbi);

    % Calculate the deformation/position of substrate, xst
    vf(n) = vu*(1-(fadh(n-1)-kmem*r(n-1))/Fm);                        % Current retrograde flow speed
    ucb = Pbi.*(vf(n)*dt+uc);
    ucsum = sum(ucb);
    us = kc*ucsum/(ks+kc*Nb);
    % Calculate the clutch force and retrograde flow
    uc = not(Pbi).*us +ucb;
    fc = kc*(uc-us);
    % xci is clutch displacement: for unbound clutches xci=xst; for bounded
    % clutches, current vf determines the next xci position
    vp(n) = vu*exp(-r(n-1)/rp);
    % vp(n) = vp0;
    r(n) = (vp(n)-vf(n))*dt+r(n-1);                                                       % Current spreading radius
    t = t+dt;
    n = n+1;

end
vp(end) = [];
r(end) = [];
vf(end) = [];

results.T = T;
results.fadh = fadh;
results.fmeant = fmeant;
results.Nbt = Nbt;
results.vp = vp;
results.vf = vf;
results.r = r;
results.Tfinal = Tfinal;
results.fmyo = (1-vf/vu)*Fm;
end