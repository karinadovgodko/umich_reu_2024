getMonskyElts = (n) -> (
    GFm = GF(2^n, Variable => a);
    returnList = {{0,1}};
    for i to 2^n-2 do (
        x = a^i;
        k = 1;
        tempsum = x^k;
        while tempsum != 0 do (
            k = k+1;
            tempsum = tempsum + x^(2^k);
        );
        returnList = append(returnList, {x, k});
    );
    return returnList;
);

getMonskyAlphas = (n, i) -> (
    myList = getMonskyElts (n);
    for elt in myList do (
        if elt#1 == i return elt#0;
    );
    return {};
);
