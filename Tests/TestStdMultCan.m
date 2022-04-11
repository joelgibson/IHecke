// There is a special Standard x Canonical -> Canonical multiplication, which is faster
// than always changing to the standard basis to multiply. Perform some double-checks
// that this faster multiplication is (A) used, and (B) returns the correct results.

SetQuitOnError(true);
SetColumns(0);
SetAssertions(3);
AttachSpec("IHecke.spec");

for cartanName in ["B2", "C3"] do
    W := CoxeterGroup(GrpFPCox, cartanName);
    HAlg := IHeckeAlgebra(W);
    H := StandardBasis(HAlg);
    C := CanonicalBasis(HAlg);

    // Test the multiplication on all pairs of elements in the group.
    Welts := EnumerateCoxeterGroup(W);
    for x in Welts, y in Welts do
        slowProduct := H.x * (H ! C.y);
        assert Parent(slowProduct) eq H;

        fastProduct := H.x * C.y;
        assert Parent(fastProduct) eq C;

        assert slowProduct eq fastProduct;
    end for;
end for;


print "TestStdMultCan passed.";
quit;
