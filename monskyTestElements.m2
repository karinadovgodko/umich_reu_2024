getMonskyElts = (m) -> (
    GFm = GF(2^m, Variable => a);
    returnList = {(0,1)};
    for i to 2^m-1 do (
        x = a^i;
        n = 2;
        tempsum = x;
        while tempsum != 0 do (
            tempsum = tempsum + x^(2^n);
            n = n+1;
        );
        returnList = append(returnList, (x, n));
    );
    return returnList;
);