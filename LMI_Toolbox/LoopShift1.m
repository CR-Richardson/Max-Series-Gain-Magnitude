% Expressing a Lurie system with ReLU nonlinearity in terms of an equivalent Lurie
% system with magnitude nonlinearity.

function new_syst = LoopShift1(syst, alpha)
A = syst.a + 0.5*alpha*syst.b*syst.c;
B = 0.5*alpha*syst.b;
C = syst.c;
D = alpha*syst.d;
new_syst = ss(A,B,C,D);
end
