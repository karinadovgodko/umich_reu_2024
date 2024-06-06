load "monskyTestElements.m2";

m=2;
alpha = getMonskyAlpha(m)
R = ring alpha [x,y,z]
g = alpha*x^2*y^2+z*(x^3+y^3+z^3+x*y*z)
r=1
n= m+r

halfq = 2^(n-1)
q=2*halfq
S = R/ ideal(x^q,y^q, z^q)
T = module S

phi = map (T,T, sub(g, S))
kerPhi = ker phi
gens kerPhi
(gens kerPhi)_(0,0)
entries basis(3*halfq-4-2^(r-1), kerPhi)
gens source basis (3*halfq-4-2^(r-1), kerPhi)
basis(3*halfq-4-2^(r-1), kerPhi)

for j from   3*halfq - 4 - 2^(r-1) to 3*halfq - 4 do (
        -- jth degree component
        returnNum = numgens source basis (j, kerPhi);
        print(returnNum);
        );
