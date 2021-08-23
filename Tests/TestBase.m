// Test creation of a Hecke algebra and its standard basis, and basic operations within the standard
// basis.

SetQuitOnError(true);
SetColumns(0);
SetAssertions(3);
AttachSpec("IHecke.spec");

// Hecke algebra creation
W := CoxeterGroup(GrpFPCox, "A2");
HAlg := IHeckeAlgebra(W);

// Hecke algebra accessors
assert Sprint(HAlg) eq "Iwahori-Hecke algebra of type A2";
assert CoxeterGroup(HAlg) eq W;
LPoly<v> := BaseRing(HAlg);

// Standard basis creation
H := IHeckeAlgebraStd(HAlg);

// Standard basis accessors
assert Sprint(H) eq "Standard basis of Iwahori-Hecke algebra of type A2, symbol H";
assert FreeModule(H) eq HAlg;
assert CoxeterGroup(H) eq W;
assert BasisSymbol(H) eq "H";
assert BasisName(H) eq "Standard basis";
assert BaseRing(H) eq BaseRing(HAlg);
assert H eq H;

// Standard basis element creation.
assert H.0 eq H.(W.0);
assert H.1 eq H.(W.1);
assert H.2 eq H.(W.2);
assert H.[1,2,1] eq H.(W ! [1,2,1]);

// Standard basis element accessors.
assert Parent(H.0) eq H;
assert Coefficient(H.0, W.0) eq 1;
assert Coefficient(H.0, W.1) eq 0;

// Standard basis coercion from scalars.
assert IsZero(H ! 0);
assert not IsZero(H ! 1);
assert H ! 2 eq 2 * H.0;

// Standard basis coercion from self.
assert H ! H.[1,2] eq H.[1,2];

// Standard basis support and arithmetic.
assert -(-H.1) eq H.1;
assert Support(H.1) eq {@ W.1 @};
assert Support(H.1 - H.1) eq {@ @};
supp, coeffs := Support(H.1 - v*H.0);
assert Setseq(supp) eq [W.0, W.1];
assert coeffs eq [-v, LPoly!1];

// Standard basis multiplication: identity.
for w in EnumerateCoxeterGroup(W) do
    assert H.0 * H.w eq H.w;
    assert H.w * H.0 eq H.w;
end for;

// Standard basis multiplication: quadratic relation.
for s in [1..Rank(W)] do
    assert H.s * H.s eq (v^-1 - v)*H.s + H.0;
end for;

// Standard basis multiplication: longest word.
assert &*[H.1, H.2, H.1] eq H.[1,2,1];

// Standard basis: formatting.
assert Sprint(H.0) eq "(1)H(id)";
assert Sprint(H ! 0) eq "(0)H(id)";
assert Sprint(H.[1,2,1]) eq "(1)H(121)";
//Sprint(H.[1,2,1] * H.1);
//assert Sprint(H.[1,2,1] * H.1) eq "(v^-1 - v)H(121) + (1)H(12)";

print "TestBase passed.";
quit;