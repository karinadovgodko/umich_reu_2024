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