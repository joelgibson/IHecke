if assigned batch then SetQuitOnError(true); else SetDebugOnError(true); end if;
SetColumns(0);
SetAssertions(3);
AttachSpec("IHecke.spec");

W := CoxeterGroup(GrpFPCox, "A2");
HAlg := IHeckeAlgebra(W);
LPoly<v> := BaseRing(HAlg);
H := IHeckeAlgebraStd(HAlg);

// Standard basis of antispherical module creation.
aH := IHeckeRModuleStd(HAlg, [2]);

// Standard basis accessors
assert Sprint(aH) eq "Standard basis of right module over Hecke algebra for Coxeter group of type A2, symbol aH, parabolic quotient [ 2 ]";
assert Parent(aH) eq HAlg;
assert CoxeterGroup(aH) eq W;
assert BasisSymbol(aH) eq "aH";
assert BasisName(aH) eq "Standard basis";
assert BaseRing(aH) eq BaseRing(HAlg);
assert Parabolic(aH) eq [2];
assert aH eq aH;

// Standard basis element creation.
assert aH.0 eq aH.(W.0);
assert aH.1 eq aH.(W.1);
assert aH.[1,2] eq aH.(W.1 * W.2);

// Action of the standard basis of the Hecke algebra on the generator of the antispherical module.
assert aH.0 * H.0 eq aH.0;
assert aH.0 * H.1 eq aH.1;
assert aH.0 * H.2 eq (-v) * aH.0;
assert aH.0 * H.[1,2] eq aH.[1,2];
assert aH.0 * H.[2,1] eq (-v) * aH.1;
assert aH.0 * H.[1,2,1] eq (-v) * aH.[1,2];


// Canonical basis of antispherical module creation.
aC := IHeckeRModuleCan(HAlg, [2]);

// Canonical basis accessors.
assert Sprint(aC) eq "Canonical basis of right module over Hecke algebra for Coxeter group of type A2, symbol aC, parabolic quotient [ 2 ]";
assert Parent(aC) eq HAlg;
assert CoxeterGroup(aC) eq W;
assert BasisSymbol(aC) eq "aC";
assert BasisName(aC) eq "Canonical basis";
assert BaseRing(aC) eq BaseRing(HAlg);
assert Parabolic(aC) eq [2];
assert aC eq aC;

// Canonical basis element creation.
assert aC.0 eq aC.(W.0);
assert aC.1 eq aC.(W.1);
assert aC.[1,2] eq aC.(W.1 * W.2);

// Action of the standard basis of the Hecke algebra on the generator of the antispherical module.
// This is the same as the quotient map.
assert aC.0 * H.0 eq aC.0;
assert aC.0 * H.1 eq aC.1;
assert aC.0 * H.2 eq 0;
assert aC.0 * H.[1,2] eq aH.[1,2];
assert aC.0 * H.[2,1] eq 0;
assert aC.0 * H.[1,2,1] eq 0;


print "TestRModBase passed.";
quit;