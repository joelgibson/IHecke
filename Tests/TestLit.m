// Test the standard functions for the Literal basis.

SetQuitOnError(true);
SetColumns(0);
SetAssertions(3);
AttachSpec("IHecke.spec");

W := CoxeterGroup(GrpFPCox, "C2");
HAlg := IHeckeAlgebra(W);
LPoly<v> := BaseRing(HAlg);
H := StandardBasis(HAlg);
C := CanonicalBasis(HAlg);
p2C := LiteralBasisInCan(HAlg, "p2C", "2-canonical basis");

// Canonical basis accessors
assert Sprint(p2C) eq "2-canonical basis of Iwahori-Hecke algebra of type B2, symbol p2C";
assert FreeModule(p2C) eq HAlg;
assert CoxeterGroup(p2C) eq W;
assert BasisSymbol(p2C) eq "p2C";
assert BasisName(p2C) eq "2-canonical basis";
assert BaseRing(p2C) eq BaseRing(HAlg);
assert p2C eq p2C;

// Set a partial canonical basis (as if we're calculating it using some algorithm).
SetBasisElement(p2C, W.0, C.0);
assert C ! p2C.0 eq C.0;

SetBasisElement(p2C, W.1, C.1);
SetBasisElement(p2C, W.2, C.2);
SetBasisElement(p2C, W![1,2], C.[1,2]);
SetBasisElement(p2C, W![2,1], C.[2,1]);
SetBasisElement(p2C, W![1,2,1], C.[1,2,1] + C.1);

assert C ! p2C.[1,2,1] eq C.[1,2,1] + C.1;
assert p2C ! C.[1,2,1] eq p2C.[1,2,1] - p2C.1;

// Check that conversion into the standard basis works.
assert H ! p2C.0 eq H.0;
assert H ! p2C.1 eq H ! C.1;
assert H ! p2C.2 eq H ! C.2;
assert H ! p2C.[1,2,1] eq H ! C.[1,2,1] + C.1;

// Check conversion out of the standard basis works.
assert p2C ! (H ! C.[1,2,1] + C.1) eq p2C.[1,2,1];

// Check serialisation.
serialised := SerialiseBasis(p2C);
text := Sprint(serialised, "Magma");
evalled := eval text;
L := DeserialiseBasis(HAlg, evalled);

// TODO: When putting (H ! ) in front of these, gets a segmentation fault.
for w in EnumerateCoxeterGroup(W) do
    if IsDefined(L, w) then
        assert H ! p2C.w eq L.w;
    end if;
end for;

// TODO: Un-comment this when it no longer crashes Magma.
// fname := GetTempDir() cat "/" cat Tempname("IHecke-TestLit-XXX");
// fd := Open(fname, "wb");
// WriteObject(fd, r);
// delete fd;
// System("rm " cat fname);


// Check unitriangular with non-identity units on diagonal works.
W := CoxeterGroup(GrpFPCox, "A1 A1");
HAlg, H, C := ShortcutIHeckeAlgebra(W);
LPoly := BaseRing(HAlg);
T := CreateLiteralBasis(HAlg, "Standard", "T", "Standard-KL basis");

for w in EnumerateCoxeterGroup(W) do
    SetBasisElement(T, w, v^(-#w) * H.w);
end for;

for w in EnumerateCoxeterGroup(W) do
    assert H ! T.w eq v^(-#w) * H.w;
    assert T ! H.w eq v^#w * T.w;
end for;

// Check that serialise/deserialise works even when we have Laurent polynomial coefficients and have renamed the
// Laurent series variable.
text := Sprint(SerialiseBasis(T), "Magma");
v := LPoly.1;
T2 := DeserialiseBasis(HAlg, eval text);

print "TestLit passed.";
quit;
