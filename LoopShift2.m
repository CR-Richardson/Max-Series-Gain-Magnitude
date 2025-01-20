% Expressing a Lurie system with magnitude nonlinearity in terms of an equivalent Lurie
% system with ReLU nonlinearity.

function new_syst = LoopShift2(syst, alpha)
A_hat = syst.a - alpha*syst.b*syst.c;
B_hat = 2*alpha*syst.b;
C_hat = syst.c;
D_hat = alpha*syst.d;
new_syst = ss(A_hat, B_hat, C_hat, D_hat);
end