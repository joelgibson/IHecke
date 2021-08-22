import "Base.m": _IHkeFModInit;

// The base ring we work over for the Hecke algebra is the Laurent polynomials.
// In Magma, a built-in structure that does this is the Laurent series.
_LaurentPolyRing<v> := LaurentSeriesRing(Integers());

////////////////////
// The Hecke algebra

declare type AlgIHke: IHkeFMod;

intrinsic IHeckeAlgebra(W::GrpFPCox) -> AlgIHke
{Create a new Hecke algebra.}
    alg := New(AlgIHke);
    name := Sprintf("Iwahori-Hecke algebra of type %o", CartanName(W));
    _IHkeFModInit(~alg, _LaurentPolyRing, name, W);
    return alg;
end intrinsic;

intrinsic 'eq'(algA::AlgIHke, algB::AlgIHke) -> BoolElt
{Hecke algebras compare equal if their Coxeter matrices (and hence groups) compare equal.}
    return algA`CoxeterMatrix eq algB`CoxeterMatrix;
end intrinsic;

intrinsic DefaultBasis(alg::AlgIHke) -> AlgIHkeStd
{}
    return IHeckeAlgebraStd(alg);
end intrinsic;


////////////////////
// The right antispherical module

declare type ASphIHke: IHkeFMod;
declare attributes ASphIHke:
    Para;   // Subsequence of [1..Rank(W)], the generators of the parabolic subgroup.

// TODO: Perhaps this should accept W as an argument rather than HAlg? Check if we need to cache anything.
intrinsic IHeckeAntiSpherical(HAlg::AlgIHke, I::SeqEnum[RngIntElt]) -> ASphIHke
{Create the right antispherical module, with basis ^I W (minimal coset representatives for the right
 cosets of W_I in W.}
    // Normalise I and ensure that it is within bounds.
    W := CoxeterGroup(HAlg);
    I := Sort(Setseq(Seqset(I)));
    require I subset [1..Rank(W)]:
        "Parabolic subset", I, "should be a subset of", [1..Rank(W)];

    asmod := New(ASphIHke);
    name := Sprintf("Antispherical module of type %o, parabolic %o", CartanName(W), I);
    _IHkeFModInit(~asmod, BaseRing(HAlg), name, W);
    asmod`Para := I;
    return asmod;
end intrinsic;

intrinsic 'eq'(asmodA::ASphIHke, asmodB::ASphIHke) -> BoolElt
{}
    return asmodA`CoxeterMatrix eq asmodB`CoxeterMatrix and asmodA`Para eq asmodB`Para;
end intrinsic;

intrinsic DefaultBasis(asmod::ASphIHke) -> ASModIHkeStd
{}
    return IHeckeAntiSphericalStd(asmod);
end intrinsic;

intrinsic Parabolic(asmod::ASphIHke) -> SeqEnum[RngIntElt]
{The parabolic generators.}
    return asmod`Para;
end intrinsic;


/////////////////////////////
// The right spherical module

declare type SphIHke: IHkeFMod;
declare attributes SphIHke:
    Para;   // Subsequence of [1..Rank(W)], the generators of the parabolic subgroup.

// TODO: Perhaps this should accept W as an argument rather than HAlg? Check if we need to cache anything.
intrinsic IHeckeSpherical(HAlg::AlgIHke, I::SeqEnum[RngIntElt]) -> SphIHke
{Create the right spherical module, with basis ^I W (minimal coset representatives for the right
 cosets of W_I in W.}
    // Normalise I and ensure that it is within bounds.
    W := CoxeterGroup(HAlg);
    I := Sort(Setseq(Seqset(I)));
    require I subset [1..Rank(W)]:
        "Parabolic subset", I, "should be a subset of", [1..Rank(W)];

    smod := New(SphIHke);
    name := Sprintf("Spherical module of type %o, parabolic %o", CartanName(W), I);
    _IHkeFModInit(~smod, BaseRing(HAlg), name, W);
    smod`Para := I;
    return smod;
end intrinsic;

intrinsic 'eq'(smodA::SphIHke, smodB::SphIHke) -> BoolElt
{}
    return smodA`CoxeterMatrix eq smodB`CoxeterMatrix and smodA`Para eq smodB`Para;
end intrinsic;

intrinsic DefaultBasis(smod::SphIHke) -> SModIHkeStd
{}
    return IHeckeSphericalStd(smod);
end intrinsic;

intrinsic Parabolic(smod::SphIHke) -> SeqEnum[RngIntElt]
{The parabolic generators.}
    return smod`Para;
end intrinsic;