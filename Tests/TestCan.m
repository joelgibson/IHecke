// Test creation of the canonical basis, and conversion to and from the standard basis.
// Also use the bar involution on the standard basis to double-check the canonical basis.

SetQuitOnError(true);
SetColumns(0);
SetAssertions(3);
AttachSpec("IHecke.spec");

// Create a Hecke algebra with a Coxeter group with an interesting canonical basis.
W := CoxeterGroup(GrpFPCox, "B3");
HAlg := IHeckeAlgebra(W);
LPoly<v> := BaseRing(HAlg);
H := IHeckeAlgebraStd(HAlg);
C := IHeckeAlgebraCan(HAlg);

// Canonical basis accessors
assert Sprint(C) eq "Canonical basis of Hecke algebra for Coxeter group of type B3, symbol C";
assert Parent(C) eq HAlg;
assert CoxeterGroup(C) eq W;
assert BasisSymbol(C) eq "C";
assert BasisName(C) eq "Canonical basis";
assert BaseRing(C) eq BaseRing(HAlg);
assert C eq C;

// Canonical basis for the type A2 parabolic dihedral subgroup on [1, 2].
assert H ! C.0 eq H.0;
assert H ! C.1 eq H.1 + v*H.0;
assert H ! C.2 eq H.2 + v*H.0;
assert H ! C.[1,2] eq H.[1, 2] + v*H.1 + v*H.2 + v^2*H.0;
assert H ! C.[2,1] eq H.[2, 1] + v*H.1 + v*H.2 + v^2*H.0;
assert H ! C.[1,2,1] eq H.[1,2,1] + v*H.[1,2] + v*H.[2,1] + v^2*H.1 + v^2*H.2 + v^3*H.0;

// Round-trip every canonical basis element through two basis conversions, and likewise for
// standard basis elements.
for w in EnumerateCoxeterGroup(W) do
    assert C ! (H ! C.w) eq C.w;
    assert H ! (C ! H.w) eq H.w;
end for;

// Bar involution on canonical basis. This fixes each basis element, and sends v -> v^-1.
for w in EnumerateCoxeterGroup(W) do
    assert Bar(C.w) eq C.w;
    assert Bar((v^2 - v^-1)*C.w) eq (v^-2 - v)*C.w;
end for;

// Bar involution on the standard basis: identity and generators.
assert Bar(H.0) eq H.0;
assert Bar(v * H.0) eq v^-1 * H.0;
assert Bar(H.1) eq H.1 + (v - v^-1)*H.0;
assert Bar(H.2) eq H.2 + (v - v^-1)*H.0;

// Check the self-duality of the whole canonical basis, using the bar involution in the standard
// basis. Bar-involution of complicated elements in the standard basis is quite slow, so we don't
// check larger groups in these short tests.
for w in EnumerateCoxeterGroup(W) do
    assert Bar(H ! C.w) eq H ! C.w;
end for;

// Check that arithmetic and equality work on mixed-basis elements.
assert H.1 + v*H.0 eq C.1;
assert H.1 + v*H.0 ne C.0;

assert C.1 - v*H.0 eq H.1;
assert IsZero(C.1 - v*H.0 - H.1);


print "TestCan passed.";
quit;