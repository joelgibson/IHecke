////////////////////
// The Hecke algebra

declare type IHkeAlg: FModIHke;

intrinsic IHeckeAlgebra(W::GrpFPCox) -> IHkeAlg
{Create a new Hecke algebra.}
    // Magma has no Laurent polynomials, but we can use the polynomial part of the series ring.
    LPolyRing<v> := LaurentSeriesRing(Integers());
    HAlg := New(IHkeAlg);
    name := Sprintf("Iwahori-Hecke algebra of type %o", CartanName(W));
    _FModIHkeInit(~HAlg, LPolyRing, name, W);
    return HAlg;
end intrinsic;

intrinsic 'eq'(HAlg1::IHkeAlg, HAlg2::IHkeAlg) -> BoolElt
{Hecke algebras compare equal if their Coxeter matrices (and hence groups) compare equal.}
    return HAlg2`CoxeterMatrix eq HAlg2`CoxeterMatrix;
end intrinsic;

intrinsic DefaultBasis(HAlg::IHkeAlg) -> IHkeAlgStd
{}
    return IHeckeAlgebraStd(HAlg);
end intrinsic;


/////////////////////////////////
// The right antispherical module

declare type IHkeASMod: FModIHke;
declare attributes IHkeASMod:
    Para;   // Subsequence of [1..Rank(W)], the generators of the parabolic subgroup.

// TODO: Perhaps this should accept W as an argument rather than HAlg? Check if we need to cache anything.
intrinsic IHeckeAntiSpherical(HAlg::IHkeAlg, I::SeqEnum[RngIntElt]) -> IHkeASMod
{Create the right antispherical module, with basis ^I W (minimal coset representatives for the right
 cosets of W_I in W.}
    // Normalise I and ensure that it is within bounds.
    W := CoxeterGroup(HAlg);
    I := Sort(Setseq(Seqset(I)));
    require I subset [1..Rank(W)]:
        "Parabolic subset", I, "should be a subset of", [1..Rank(W)];

    asmod := New(IHkeASMod);
    name := Sprintf("Antispherical module of type %o, parabolic %o", CartanName(W), I);
    _FModIHkeInit(~asmod, BaseRing(HAlg), name, W);
    asmod`Para := I;
    return asmod;
end intrinsic;

intrinsic 'eq'(ASMod1::IHkeASMod, ASMod2::IHkeASMod) -> BoolElt
{}
    return ASMod1`CoxeterMatrix eq ASMod2`CoxeterMatrix and ASMod1`Para eq ASMod2`Para;
end intrinsic;

intrinsic DefaultBasis(ASMod::IHkeASMod) -> ASModIHkeStd
{}
    return IHeckeAntiSphericalStd(ASMod);
end intrinsic;

intrinsic Parabolic(ASMod::IHkeASMod) -> SeqEnum[RngIntElt]
{The parabolic generators.}
    return ASMod`Para;
end intrinsic;


/////////////////////////////
// The right spherical module

declare type IHkeSMod: FModIHke;
declare attributes IHkeSMod:
    Para;   // Subsequence of [1..Rank(W)], the generators of the parabolic subgroup.

// TODO: Perhaps this should accept W as an argument rather than HAlg? Check if we need to cache anything.
intrinsic IHeckeSpherical(HAlg::IHkeAlg, I::SeqEnum[RngIntElt]) -> IHkeSMod
{Create the right spherical module, with basis ^I W (minimal coset representatives for the right
 cosets of W_I in W.}
    // Normalise I and ensure that it is within bounds.
    W := CoxeterGroup(HAlg);
    I := Sort(Setseq(Seqset(I)));
    require I subset [1..Rank(W)]:
        "Parabolic subset", I, "should be a subset of", [1..Rank(W)];

    smod := New(IHkeSMod);
    name := Sprintf("Spherical module of type %o, parabolic %o", CartanName(W), I);
    _FModIHkeInit(~smod, BaseRing(HAlg), name, W);
    smod`Para := I;
    return smod;
end intrinsic;

intrinsic 'eq'(SMod1::IHkeSMod, SMod2::IHkeSMod) -> BoolElt
{}
    return SMod1`CoxeterMatrix eq SMod2`CoxeterMatrix and SMod1`Para eq SMod2`Para;
end intrinsic;

intrinsic DefaultBasis(SMod::IHkeSMod) -> SModIHkeStd
{}
    return IHeckeSphericalStd(SMod);
end intrinsic;

intrinsic Parabolic(SMod::IHkeSMod) -> SeqEnum[RngIntElt]
{The parabolic generators.}
    return SMod`Para;
end intrinsic;