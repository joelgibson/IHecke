if assigned batch then SetQuitOnError(true); else SetDebugOnError(true); end if;
SetColumns(0);
SetAssertions(3);
AttachSpec("IHecke.spec");

// Run an S3 example
W := CoxeterGroup(GrpFPCox, "A2");
HAlg := IHeckeAlgebra(W);
C := IHeckeAlgebraCan(HAlg);

// Enumerating Coxeter groups
assert #EnumerateCoxeterGroup(W: lengthBound := 0) eq 1;
assert #EnumerateCoxeterGroup(W: lengthBound := 1) eq 3;
assert #EnumerateCoxeterGroup(W: lengthBound := 2) eq 5;
assert #EnumerateCoxeterGroup(W: lengthBound := 3) eq 6;
assert #EnumerateCoxeterGroup(W: lengthBound := 4) eq 6;
assert #EnumerateCoxeterGroup(W) eq 6;

// Cell computations: in A2 we have 3 two-sided cells, and 4 left or right-sided cells.
left, right, two := Cells(C);
assert #Components(left) eq 4;
assert #Components(right) eq 4;
assert #Components(two) eq 3;

// Test cell order
assert CellLe(left, W.0, W.1);
assert CellLe(left, W.0, W.2);
assert CellLe(left, W.0, LongestElement(W));
assert not CellLe(left, W.1, W.0);
assert not CellLe(left, W.2, W.0);


print "TestCells passed.";
quit;