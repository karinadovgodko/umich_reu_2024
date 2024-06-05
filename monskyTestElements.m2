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
getMonskyAlpha = (n, i) -> (
    if not (instance(i, ZZ) and i>0 and instance(n, ZZ) and n>0 and i<= 2*n) then error "Inputs must be positive integers with the second at most twice the first";
    privateGF2n = GF(2^n, Variable => a);
    for j to 2^n-2 do (

        elt = a^j;
        k = 0;
        tempsum = elt;

        while tempsum != 0 do (
            k += 1;
            tempsum += elt^(2^k);
        );

        if k+1 == i then return elt;

    );
    return sub(0, privateGF2n);
);
