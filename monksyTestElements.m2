getMonskyElts = (n) -> (
    GFn = GF(2^n, Variable => a);
    returnList = {{0,1}};
    for i to 2^n-2 do (
        elt = a^i;
        k = 0;
        tempsum = elt;


        while tempsum != 0 do (
            k = k+1;
            --print "(" |x|"," | tempsum |")"; 
            -- print (x, tempsum, tempsum != 0);
            tempsum += elt^(2^k);
            if k == 50 then break;
        );
    --    print (x, tempsum, tempsum != 0);
        returnList = append(returnList, {elt, k+ 1});

    );
    return returnList;
);

getMonskyAlphas = (n, i) -> (
    myList = getMonskyElts (n);
    for elt in myList do (
        if elt#1 == i then return elt#0;
    );
    return myList#0#0;
);
