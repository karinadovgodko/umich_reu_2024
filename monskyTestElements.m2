-- since the order of the multiplicative group of GF(2,n) is 2^n-1, an element has finite m if and only if it does for m up to 2n.
getMonskyElts = (n) -> (
    privateGF2n = GF(2^n, Variable => a);
    returnList = {{sub(0,privateGF2n),1}};
    for i to 2^n-2 do (
        elt = a^i;
        k = 0;
        tempsum = elt;


        while tempsum != 0 do (
            k += 1;
            tempsum += elt^(2^k);
        );
        returnList = append(returnList, {elt, k+ 1});

    );
    return returnList;
);

-- note: This method will return 0 in the field if there is no such element because we need something comparable to a general element of the ring and we are generally not interested in 0. Also m\le 2n always
getMonskyAlpha = (m) -> (
    -- assuming there is an element for any m in some field with 2^n elements
    if not (instance(m, ZZ) and m>0 ) then error "Inputs must be positive integers with the second at most twice the first";

    if m == 1 then return sub(0, GF(2)); 
    n = 1; 
    while true do (
        privateGF2n = GF(2^n, Variable => a);
        for j to 2^n-2 do (

            elt = a^j;
            k = 0;
            tempsum = elt;

            while tempsum != 0 do (
                k += 1;
                tempsum += elt^(2^k);
            );

            if k+1 == m then return elt;

        );
        n += 1; 
    )
);

getTwoMonskyElementsInSameField = (m1, m2) -> (

n = 1;
while true do (
    alpha = 0;
    beta = 0;

    listPairs = getMonskyElts(n);

    for pair in listPairs do (
        if pair#1 == m1 then alpha = pair#0; 
        if pair#1 == m2 then beta = pair#0;
    );

    if alpha != 0 and beta != 0 then return {alpha, beta};
    n += 1; 
);
);