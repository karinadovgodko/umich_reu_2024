hkmApproximation = (f,n) -> (
    if not (instance(f, RingElement))then error "Input1 must be a polynomial";
    if not (instance (n, ZZ) and n>0) then error "Input2 must be a positive integer";
    S = ring f;
    Schar = char(S);

    if not (isPrime(Schar)) then error "Input1 must have coefficients in prime characteristic";

    Sgens = gens S;
    powernumber = Schar^n;
    myList = {f};
     for term in Sgens do (
        myList = append(myList, term^powernumber);
    );
    myIdeal = ideal(myList);
    R = module(S)/module(myIdeal);


    --Note: I believe the a module is finite dimensional iff it has finite length. I don't think this is a necessary error catch because R should always be finite-dimensional

    --if not isFinite length R then error "The quotient of ring input1 by the ideal is not a finite dimensional vector space";

    return (degree R)/(Schar^(n*(numgens S-1)));
);

hkmApproximationComponent = (f, n, i) -> (
    if not (instance(f, RingElement) )then error "Input1 must be a polynomial";
    if not (instance (n, ZZ) and n>0) then error "Input2 must be a positive integer";
    S = ring f;
    Schar = char(S);

    if not (isPrime(Schar)) then error "Input1 must have coefficients in prime characteristic";

    Sgens = gens S;
    powernumber = Schar^n;
    myList = {f};
     for term in Sgens do (
        myList = append(myList, term^powernumber);
    );
    myIdeal = ideal(myList);
    R = module(S)/module(myIdeal);
    


    --Note: I believe the a module is finite dimensional iff it has finite length. I don't think this is a necessary error catch because R should always be finite-dimensional

    --if not isFinite length R then error "The quotient of ring input1 by the ideal is not a finite dimensional vector space";

    --also, for some reason numcols basis (i,R) is SIGNIFICANTLY faster than hilbertFunction(i,R). However, calling the hkmApproximationComponent iteratively is slower than the main method.

    return numcols basis (i, R);
);

algebraicHKExample = (f, n) -> (
    if not (instance(f, RingElement) and (numgens ring f) == 1) or ideal((gens ring f)#0) == ring f then error "Input1 must be a polynomial in one variable";
    if (degree(f) == {}) then error "Input1 cannot be constant";
    if not (instance(n, ZZ) and n>0) then error "Input2 must be a positive integer";
    if ((degree(f))#0 == 1 and (f(1) == 0 or f(0) == 0 or f(-1) == 0)) then error "Input1 cannot be of the form x-1, x, or x+1";

    S = ring f;
    gensS = gens S;
    charS = char S;
    factors = factor f;

    if not (factors#0#1 == 1 and (#factors == 1 or (#factors == 2 and ring f == ideal(factors#1#0)))) then error "Input1 must be irreducible";
    if not (isPrime charS and isField coefficientRing S) then error "Input1 must be a polynomial over a field of prime characteristic";

    k = S / ideal(f);
    alpha = sub (gensS#0, k);
    K = k[X,Y];
    qPoly = Y*(Y^2-X^2)*(Y-alpha*X);
    return hkmApproximation(qPoly, n);

)

transcendentalHKExample = (R, n) -> (
    if not (instance(n, ZZ) and n>0) then error "Input2 must be a positive integer";
    if not (isField R and isPrime(char R)) then error "Input1 must be a field of positive characteristic";
    S = R[t];
    k = frac S;
    K = k[X,Y];
    qPoly = Y*(Y^2-X^2)*(Y-t*X);
    return hkmApproximation(qPoly, n);
)