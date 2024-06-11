load "monskyTestElements.m2"


baseSegreRing = (f, g) -> (
    if not (instance(f, RingElement) and instance(g, RingElement) and isHomogeneous(f) and isHomogeneous(g)) then error "Inputs must be homogeneous polynomials"; 

    R = ring f;
    S = ring g;
    Q1 = R / ideal(f);
    Q2 = S / ideal(g);

    if not coefficientRing R === coefficientRing S then error "Polynomials must have coefficients in the same ring";
    PolynomialSegreRing = coefficientRing R[z_(1,1) .. z_(numgens Q1, numgens Q2)];
    return PolynomialSegreRing;
);

segreProduct = (f, g) -> (
    -- Ensure the inputs are valid
    if not (instance(f, RingElement) and instance(g, RingElement) and isHomogeneous(f) and isHomogeneous(g)) then error "Inputs must be homogeneous polynomials"; 

    R = ring f;
    S = ring g;

    if not coefficientRing R === coefficientRing S then error "Polynomials must have coefficients in the same ring";

    Q1 = R / ideal(f);
    Q2 = S / ideal(g);
    SegreRing = baseSegreRing(f, g);

    -- Compute the tensor product ring
    tensorRing = Q1 ** Q2;
    tensorGens = gens tensorRing;
    SegreGens = gens SegreRing;


    -- Define the mapping for the ring homomorphism
    relationsList = {};
    for i to numgens SegreRing - 1 do (
        numberuno = i // numgens Q2;
        numberdos = numgens Q1 + (i % numgens Q2);
        relationsList = append( relationsList, tensorGens#numberuno*tensorGens#numberdos);
    );


    -- Define the ring homomorphism from the tensor product to the Segre product ring
    relationsMap = map(tensorRing, SegreRing, relationsList);

    -- Compute the kernel of the homomorphism
    kerRelationsMap = kernel relationsMap;

    -- Define the Segre product as the quotient of the tensor product by this kernel
    outputSegreProductPresentation = SegreRing / kerRelationsMap;

    return outputSegreProductPresentation;
);


segreProductMultiplicity(m1, m2, n) -> (
    alphas = getTwoMonskyElementsInSameField(m1, m2);
    alpha = alphas#0; 
    beta = alphas#1; 
    myRing = (ring alpha)[x_1..x_3];
    g1 = alpha*x_1^2*x_2^2+x_3*(x_1^3+x_2^3+x_3^3+x_1*x_2*x_3);
    g2 = beta*x_1^2*x_2^2+x_3*(x_1^3+x_2^3+x_3^3+x_1*x_2*x_3);
    
    segreRing = segreProduct(g1, g2);

    Sgens = gens segreRing;
    q = 2^n;
    myList = {};
     for term in Sgens do (
        myList = append(myList, term^q);
    );
    myIdeal = ideal(myList);
    R = module(segreRing)/module(myIdeal);

    d = dim segreRing; 

    return (degree R)/(2^(n*d));
);
