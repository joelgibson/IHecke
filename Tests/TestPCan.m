// Test loading a few p-canonical basis elements, and basis conversions.

SetQuitOnError(true);
SetColumns(0);
SetAssertions(3);
AttachSpec("IHecke.spec");

for cartanName in ["B2", "C3"] do
    // Loading p-canonical bases
    W := CoxeterGroup(GrpFPCox, cartanName);
    HAlg := IHeckeAlgebra(W);
    H := StandardBasis(HAlg);
    C := CanonicalBasis(HAlg);
    pC := IHeckeAlgebraPCan(HAlg, cartanName, 2: quiet := true);

    // p-canonical basis accessors
    assert FreeModule(pC) eq HAlg;
    assert CoxeterGroup(pC) eq W;
    assert BasisSymbol(pC) eq "p2C";
    assert BasisName(pC) eq "2-canonical basis";
    assert BaseRing(pC) eq BaseRing(HAlg);
    assert pC eq pC;

    // Round-trip through various basis conversions.
    for w in EnumerateCoxeterGroup(W) do
        assert pC ! (C ! pC.w) eq pC.w;
        assert pC ! (H ! pC.w) eq pC.w;

        assert H ! (pC ! H.w) eq H.w;
        assert C ! (pC ! C.w) eq C.w;
    end for;
end for;


print "TestPCan passed.";
quit;