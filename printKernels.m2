load "HKMApproximation.m2"

printKernels = (name, n1, n2) -> (
    f = name << ""; 
    for n from n1 to n2 do (
        dimensions =  mutableMatrix (ZZ, 3*2^5 - 2, 3);
f <<  "------------ GF 2^" << n  <<  "--------------" << endl;
    for i from 3 to 5  do (
    -- ith Frobenius iteration 
	for j from 0 to 3*2^i - 3 do (
    -- jth degree component
    dimensions_(j, i-3) = hkmKernelIndex(n, n, i, j)
    );
);
    f << dimensions << endl;
        );

    f << close;
    );
  