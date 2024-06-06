load "HKMApproximation.m2";

--name is the name of the file. The method constructs tables for alpha in GF(2^i) for i from n1 to n2 with m(alpha)=i. Then for j from n3 to n4, we also compute the j-th Frobenius iterates.
printKernels = (name, m1, m2, upperFrob) -> (
    -- (m1, m2) range of m(alpha) 
    if not (upperFrob > m2) then error "more frobenius iterations needed";
    f = name << ""; 
        f << "Rows are dimensions of kernel from 
        3*2^(i-1) - 4 - 2^(r-1) to 3*2^(i-1) - 4 for i the number of Frobenius iterations and r = i - m.
        Columns are Frobenius iterations" << endl;
    for m from m1 to m2 do (
        f <<  "------------ m = " << m  <<  "--------------" << endl;
    dimensions =  mutableMatrix (ZZ, 2^(upperFrob-m-1) + 1, upperFrob -  m);
    for i from m + 1 to upperFrob do (
     -- ith Frobenius iteration 
        r = i - m;
        for j from   3*2^(i-1) - 4 - 2^(r-1) to 3*2^(i-1) - 4 do (
        -- jth degree component
        dimensions_(j - (3*2^(i-1) - 4 - 2^(r-1)), i-(m+1)) = hkmKernelIndex(m, i, j);
        );
    );
    f << dimensions << endl;
        );

    f << close;
    );
  
  --name is the name of the file. The method constructs tables for alpha in GF(2^i) for i from n1 to n2 with m(alpha)=i. Then for j from n3 to n4, we also compute the j-th Frobenius iterates.
printKernelsRange = (name, m1, m2, upperFrob, range) -> (
    -- (m1, m2) range of m(alpha) 
    if not (upperFrob > m2) then error "more frobenius iterations needed";
    f = name << ""; 
        f << "Rows are dimensions of kernel from 
        3*2^(i-1) - 4 - 2^(r-1) to 3*2^(i-1) - 4 for i the number of Frobenius iterations and r = i - m.
        Columns are Frobenius iterations" << endl;
    for m from m1 to m2 do (
        f <<  "------------ m = " << m  <<  "--------------" << endl;
        tempNum = min(range+m, upperFrob);
    dimensions =  mutableMatrix (ZZ, 2^(tempNum-m-1) + 1, tempNum -  m);
    for i from m + 1 to tempNum do (
     -- ith Frobenius iteration 
        r = i - m;
        for j from   3*2^(i-1) - 4 - 2^(r-1) to 3*2^(i-1) - 4 do (
        -- jth degree component
        dimensions_(j - (3*2^(i-1) - 4 - 2^(r-1)), i-(m+1)) = hkmKernelIndex(m, i, j);
        );
    );
    f << dimensions << endl;
        );

    f << close;
    );
  
  --name is the name of the file. The method gives basis elements for the kernels of multiplication by  g_alpha with m(alpha)=m for frobenius iterates from m+1 to upperFrob.
printKernelsBasis = (name, m , upperFrob) -> (
    if not (upperFrob > m) then error "more frobenius iterations needed";
    f = name << ""; 

  
    for i from m + 1 to upperFrob do (
        -- ith Frobenius iteration 
        f <<  "------------ Frobenius Iteration = " << i  <<  "--------------" << endl;

        r = i - m;
        for j from   3*2^(i-1) - 4 - 2^(r-1) to 3*2^(i-1) - 4 do (
            
            -- jth degree component
            base = hkmKernelIndexBasis(m, i, j);

            f << base << endl;
        );
    );
    f << close;
    );
  
