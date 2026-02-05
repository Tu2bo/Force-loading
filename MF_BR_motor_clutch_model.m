function dydt = MF_BR_motor_clutch_model(t, y, kc, ks, vu, kon0, koff0, Fb, Nc, Fm, kmem, rp, vbeta, alpha)
% Mean-field approximation of the BR-integrated motor–clutch model
%
% y(1): total adhesion force
% y(2): number of bound clutches
% y(3): cell spreading radius
dydt = zeros(3,1);
vf = vu*(1-(y(1)-kmem*y(3))/Fm);
vp = vu*exp(-y(3)/rp);
koff = min(koff0*exp(alpha*y(1)/y(2)/Fb),1e6);
roff = koff/koff0;


dydt(1) = ks*y(2)*kc./(ks+y(2)*kc)*(vf-vbeta*roff*log(roff));   %dfadh/dt
dydt(2) = (Nc-y(2))*kon0-y(2)*koff;                                         %dNb/dt
dydt(3) = vp-vf;                                                                       %dr/dt
end