load "HKMApproximation.m2";

--name is the name of the file. The method constructs tables for alpha in GF(2^i) for i from n1 to n2 with m(alpha)=i. Then for j from n3 to n4, we also compute the j-th Frobenius iterates.
printKernels = (name, n1, n2, n3, n4) -> (
    f = name << ""; 
    for n from n1 to n2 do (
        dimensions =  mutableMatrix (ZZ, 3*2^n4 - 2, n4-n3+1);
f <<  "------------ GF 2^" << n  <<  "--------------" << endl;
    for i from n3 to n4  do (
        -- ith Frobenius iteration 
        for j from 0 to 3*2^i - 3 do (
        -- jth degree component
        dimensions_(j, i-n3) = hkmKernelIndex(n, n, i, j)
        );
    );
    f << dimensions << endl;
        );

    f << close;
    );
  
