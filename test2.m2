load "HKMApproximation.m2";

--load "segreProductMethod.m2";


l=2;
r=4;
n=2*l+r;
halfq=2^(n-1);
q=2*halfq;

myField = GF(2, 2, Variable => a);
alpha = a;
myPolyRing = myField [x,y,z];
Ax = x^2+y*z;
Ay = y^2+x*z;
h = alpha*z^4+Ax * Ay;

myIdeal = ideal(x^q,y^q,z^q);

R = myPolyRing/myIdeal;

--  dim (myPolyRing/ideal(h))

phi = map(R^1,R^1, h);

kerPhi = ker phi;

u = (gens kerPhi)_(0,0)

basisTestElt = (a,b,c,i) -> (
    if not b== 0 or b==1 then error "b must be 0 or 1";
    return sub((Ax+Ay)^i*x^a*y^b*z^c*u, R);
);

basisTestElt(1,0,0,0)

-- cokerPhi = coker phi;






-- mySum = 0;

-- for i to 3 * q-3 do (
--     mySum += numgens source basis (i, kerPhi);
-- )

-- mySum

-- 3*4^n+4^r
-- mySum -( 3*4^n+4^r)

-- for s to 3+2^(r-1) do (
--     print ("s = "| s , " dimension of kernel ="| numgens source basis (3*halfq-7-2^(r-1)+s, kerPhi));
-- );
