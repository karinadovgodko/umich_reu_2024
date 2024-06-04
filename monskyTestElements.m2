getMonskyElts = (m) -> (
    GFm = GF(2^m, Variable => a);
    returnList = {(0,1)};
    for i to 2^m-1 do (
        x = a^i;
        n = 1;
        tempsum = x^n;
        while tempsum != 0 do (
            n = n+1;
            tempsum = tempsum + x^(2^n);
        );
        returnList = append(returnList, (x, n));
    );
    return returnList;
);
