%Ruihao Xue 2025
%Mean-field approximation of the BR-integrated motor-clutch model
function dydT = MFSingleSpringEntDim(t, y, K, V, Vp0, Fp, Ron, Fm, Kmem, alpha)
% y(1) F
% y(2) P
% y(3) R
dydT = zeros(3,1);

y(2)= max(y(2), 1e-4);
Roff = min(exp(alpha*y(1)/y(2)), 1e6);
H = 1-(y(1)-Kmem*y(3))/Fm;


dydT(1) = K*(V*H-Roff*log(Roff));            %dF/dT
dydT(2) = (1-y(2))*Ron-y(2)*Roff;            %dP/dT
dydT(3) = Vp0*(exp(-y(3)/Fp)-H) ;           %dR/dT


end
