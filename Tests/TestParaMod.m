// Test creation of the spherical and antispherical modules, and the action of the Hecke algebra
// upon these modules.

SetQuitOnError(true);
SetColumns(0);
SetAssertions(3);
AttachSpec("IHecke.spec");

// Return true if the zero-argument function 'fn' throws an error when called.
// Used for making sure that non-I-minimal elements are rejected.
function throwsError(fn)
    try
        _ := fn();
        return false;
    catch e
        return true;
    end try;
end function;

// Hecke algebra setup.
W := CoxeterGroup(GrpFPCox, "A2");
HAlg := IHeckeAlgebra(W);
LPoly<v> := BaseRing(HAlg);
H := IHeckeAlgebraStd(HAlg);
C := IHeckeAlgebraCan(HAlg);


// Antispherical module creation.
ASMod := IHeckeAntiSpherical(HAlg, [2]);

// Module accessors
assert Sprint(ASMod) eq "Antispherical module of type A2, parabolic [ 2 ]";
assert CoxeterGroup(ASMod) eq W;

// Standard basis creation
aH := IHeckeAntiSphericalStd(ASMod);

// Standard basis accessors
assert Sprint(aH) eq "Standard basis of Antispherical module of type A2, parabolic [ 2 ], symbol aH";
assert Parent(aH) eq ASMod;
assert CoxeterGroup(aH) eq W;
assert BasisSymbol(aH) eq "aH";
assert BasisName(aH) eq "Standard basis";
assert BaseRing(aH) eq BaseRing(HAlg);
assert aH eq aH;

// Standard basis element creation.
assert aH.0 eq aH.(W.0);
assert aH.1 eq aH.(W.1);
assert throwsError(func< | aH.2 >);
assert aH.[1,2] eq aH.(W.1 * W.2);
assert throwsError(func< | aH.[2,1]>);
assert throwsError(func< | aH.[1,2,1]>);

// Should not be able to coerce scalars, aside from zero.
assert IsZero(aH ! 0);
assert throwsError(func< | aH ! 1>);

// "Tensor surjection" of the standard basis of the Hecke algebra onto the antispherical module.
assert aH.0 * H.0 eq aH.0;
assert aH.0 * H.1 eq aH.1;
assert aH.0 * H.2 eq (-v) * aH.0;
assert aH.0 * H.[1, 2] eq aH.[1, 2];
assert aH.0 * H.[2, 1] eq (-v) * aH.1;
assert aH.0 * H.[2, 1, 2] eq (-v) * aH.[1, 2];

// "Tensor surjection" of the canonical basis of the Hecke algebra onto the antispherical module.
assert aH.0 * C.0 eq aH.0;
assert aH.0 * C.1 eq aH.1 + v*aH.0;
assert IsZero(aH.0 * C.2);
assert aH.0 * C.[1, 2] eq aH.[1, 2] + v*aH.1;
assert IsZero(aH.0 * C.[2, 1]);
assert IsZero(aH.0 * C.[2, 1, 2]);


// Spherica; module creation.
SMod := IHeckeSpherical(HAlg, [2]);

// Module accessors
assert Sprint(SMod) eq "Spherical module of type A2, parabolic [ 2 ]";
assert CoxeterGroup(SMod) eq W;

// Standard basis creation
sH := IHeckeSphericalStd(SMod);

// Standard basis accessors
assert Sprint(sH) eq "Standard basis of Spherical module of type A2, parabolic [ 2 ], symbol sH";
assert Parent(sH) eq SMod;
assert CoxeterGroup(sH) eq W;
assert BasisSymbol(sH) eq "sH";
assert BasisName(sH) eq "Standard basis";
assert BaseRing(sH) eq BaseRing(HAlg);
assert sH eq sH;

// Standard basis element creation.
assert sH.0 eq sH.(W.0);
assert sH.1 eq sH.(W.1);
assert throwsError(func< | sH.2 >);
assert sH.[1,2] eq sH.(W.1 * W.2);
assert throwsError(func< | sH.[2,1]>);
assert throwsError(func< | sH.[1,2,1]>);

// Should not be able to coerce scalars, aside from zero.
assert IsZero(sH ! 0);
assert throwsError(func< | sH ! 1>);

// "Tensor surjection" of the standard basis of the Hecke algebra onto the spherical module.
assert sH.0 * H.0 eq sH.0;
assert sH.0 * H.1 eq sH.1;
assert sH.0 * H.2 eq v^-1 * sH.0;
assert sH.0 * H.[1, 2] eq sH.[1, 2];
assert sH.0 * H.[2, 1] eq v^-1 * sH.1;
assert sH.0 * H.[2, 1, 2] eq v^-1 * sH.[1, 2];

// "Tensor surjection" of the canonical basis of the Hecke algebra onto the spherical module.
assert sH.0 * C.0 eq sH.0;
assert sH.0 * C.1 eq sH.1 + v*sH.0;
assert sH.0 * C.2 eq (v + v^-1)*sH.0;
assert sH.0 * C.[1, 2] eq sH.[1, 2] + v*sH.1 + (1 + v^2)*sH.0;
assert sH.0 * C.[2, 1] eq (v^-1 + v)*sH.1 + (1 + v^2)*sH.0;
assert sH.0 * C.[2, 1, 2] eq (v + v^-1)*sH.[1,2] + (1 + v^2)*sH.1 + (v + v^3)*sH.0;


print "TestParaMod passed.";
quit;