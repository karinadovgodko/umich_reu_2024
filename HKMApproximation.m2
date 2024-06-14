load "monskyTestElements.m2";


hkmComponentCokernel = (f, n, i) -> (
    if not (instance(f, RingElement) )then error "Input1 must be a polynomial";
    if not (instance (n, ZZ) and n>0) then error "Input2 must be a positive integer";
    S = ring f;
    Schar = char(S);

    if not (isPrime(Schar)) then error "Input1 must have coefficients in prime characteristic";

    Sgens = gens S;
    powernumber = Schar^n;
    myList = {};
     for term in Sgens do (
        myList = append(myList, term^powernumber);
    );
    myIdeal = ideal(myList);
    R =S/myIdeal;

    phi = map (R^1, R^1, f);

    myCokernel = cokernel phi;
    


    --Probably not needed here: if not isFinite length R then error "The quotient of ring input1 by the ideal is not a finite dimensional vector space";

    --also, for some reason numcols basis (i,R) is SIGNIFICANTLY faster than hilbertFunction(i,R). However, calling the hkmApproximationComponent iteratively is slower than the main method.

    return numgens source basis (i, myCokernel);
);



hkmComponentKernel = (f, n, i) -> (
    -- compute dimension of ith degree component of kernel of multiplication by f (quotienting by char^n th powers of generators)
    if not (instance(f, RingElement) )then error "Input1 must be a polynomial";
    if not (instance (n, ZZ) and n>0) then error "Input2 must be a positive integer";
    S = ring f;
    Schar = char(S);

    if not (isPrime(Schar)) then error "Input1 must have coefficients in prime characteristic";

    Sgens = gens S;
    powernumber = Schar^n;
    myList = {};
     for term in Sgens do (
        myList = append(myList, term^powernumber);
    );


    myIdeal = ideal(myList);
    R = S/myIdeal;
    T = module R;

    phi = map(T, T, sub(f, R));

    phiKer = kernel phi;

    return numgens source basis(i, phiKer);
);

hkmComponentKernelBasis = (f, n, i) -> (
    -- compute dimension of ith degree component of kernel of multiplication by f (quotienting by char^n th powers of generators)
    if not (instance(f, RingElement) )then error "Input1 must be a polynomial";
    if not (instance (n, ZZ) and n>0) then error "Input2 must be a positive integer";
    S = ring f;
    Schar = char(S);

    if not (isPrime(Schar)) then error "Input1 must have coefficients in prime characteristic";

    Sgens = gens S;
    powernumber = Schar^n;
    myList = {};
     for term in Sgens do (
        myList = append(myList, term^powernumber);
    );


    myIdeal = ideal(myList);
    R = S/myIdeal;
    T = module R;

    phi = map(T, T, sub(f, R));

    phiKer = ker phi;

    return basis(i, phiKer);
);

hkmKernelAlpha = (alpha, n, i) -> (
    if not (instance (alpha, RingElement) and char (ring alpha) ==2) then error "alpha must be in some GF(2,k)";
    myRing = (ring alpha)[x,y,z];
    g = alpha*x^2*y^2+z*(x^3+y^3+z^3+x*y*z);
    return hkmComponentKernel(g, n, i);
);

hkmKernelAlphaBasis = (alpha, n, i) -> (
    if not (instance (alpha, RingElement) and char (ring alpha) ==2) then error "alpha must be in some GF(2,k)";
    myRing = (ring alpha)[x,y,z];
    g = alpha*x^2*y^2+z*(x^3+y^3+z^3+x*y*z);
    return hkmComponentKernelBasis(g, n, i);
);

hkmKernelIndex = (m, n , i) -> (
    --m as in m(alpha), quotienting by 2^nth powers of generators. returns i th degree part of kernel
    alpha = getMonskyAlpha(m);
    return hkmKernelAlpha(alpha, n, i);
);

hkmKernelIndexBasis = (m, n , i) -> (
    --m as in m(alpha), quotienting by 2^nth powers of generators. returns i th degree part of kernel
    alpha = getMonskyAlpha(m);
    return hkmKernelAlphaBasis(alpha, n, i);
);


cokerPrediction = (m, r, i) -> (
    q = 2^(m+r);
    if not (m <= r and r >= 3) then error "r must be at least as big as m, r must be at least 3";
    returnValue = 0;
    if (-4 <= i and i < 0) then returnValue = (i + 6) * (i + 5) / 2;
     if (0 <= i and i < q - 4) then returnValue = 4 * i + 14;
     if (q - 4 <= i and i < q) then returnValue = 4 * i + 14 - 3/2 * (i - q + 6) * (i - q + 5);
     if (q <= i and i < 3 * q / 2 - 4 - 2^(r - 1)) then returnValue = 12 * q - 8 * i - 28;
    if (3 * q / 2 - 4 - 2^(r - 1) <= i and i < 3 * q / 2 - 2^(r - 1)) then returnValue = 12 * q - 8 * i - 28 + (i - 3 * q / 2 + 4 + 2^(r - 1) + 2) * (i - 3 * q / 2 + 4 + 2^(r - 1) + 1) / 2;
    if (3 * q / 2 - 2^(r - 1) <= i and i < 3 * q / 2 - 6 + 2^(r - 1)) then returnValue = (6 + 2^(1 - m)) * q - 4 * i - 14;
    if (3 * q / 2 - 6 + 2^(r - 1) <= i and i < 3 * q / 2 - 2 + 2^(r - 1)) then returnValue = (3 * q / 2 - 1 + 2^(r - 1) - i) * (3 * q / 2 - 2 + 2^(r - 1) - i) / 2;
    return sub(returnValue,ZZ);
);


-- algebraicHKExample = (f, n) -> (
--     if not (instance(f, RingElement) and (numgens ring f) == 1) or ideal((gens ring f)#0) == ring f then error "Input1 must be a polynomial in one variable";
--     if (degree(f) == {}) then error "Input1 cannot be constant";
--     if not (instance(n, ZZ) and n>0) then error "Input2 must be a positive integer";
--     if ((degree(f))#0 == 1 and (f(1) == 0 or f(0) == 0 or f(-1) == 0)) then error "Input1 cannot be of the form x-1, x, or x+1";

--     S = ring f;
--     gensS = gens S;
--     charS = char S;
--     factors = factor f;

--     if not (factors#0#1 == 1 and (#factors == 1 or (#factors == 2 and ring f == ideal(factors#1#0)))) then error "Input1 must be irreducible";
--     if not (isPrime charS and isField coefficientRing S) then error "Input1 must be a polynomial over a field of prime characteristic";

--     k = S / ideal(f);
--     alpha = sub (gensS#0, k);
--     -- alpha is x subject to the relation f in k  
--     K = k[X,Y];
--     qPoly = Y*(Y^2-X^2)*(Y-alpha*X);
--     return hkmApproximation(qPoly, n);

-- )

-- transcendentalHKExample = (R, n) -> (
--     if not (instance(n, ZZ) and n>0) then error "Input2 must be a positive integer";
--     if not (isField R and isPrime(char R)) then error "Input1 must be a field of positive characteristic";
--     S = R[t];
--     k = frac S;
--     K = k[X,Y];
--     qPoly = Y*(Y^2-X^2)*(Y-t*X);
--     return hkmApproximation(qPoly, n);
-- )
