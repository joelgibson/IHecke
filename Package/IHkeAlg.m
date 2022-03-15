////////////////////
// The Hecke algebra

declare type IHkeAlg: FModIHke;

intrinsic IHeckeAlgebra(W::GrpFPCox) -> IHkeAlg
{Create a new Hecke algebra.}
    // Magma has no Laurent polynomials, but we can use the polynomial part of the series ring.
    LPolyRing<v> := LaurentSeriesRing(Integers());
    HAlg := New(IHkeAlg);
    name := Sprintf("Iwahori-Hecke algebra of type %o", _IHkeCartanName(W));
    _FModIHkeInit(~HAlg, LPolyRing, name, W);
    return HAlg;
end intrinsic;

intrinsic 'eq'(HAlg1::IHkeAlg, HAlg2::IHkeAlg) -> BoolElt
{Hecke algebras compare equal if they are instantiated from the same Coxeter group object.}
    return CoxeterGroup(HAlg1) cmpeq CoxeterGroup(HAlg2);
end intrinsic;


/////////////////////////////////
// The right antispherical module

declare type IHkeASMod: FModIHke;
declare attributes IHkeASMod:
    Para;   // Subsequence of [1..Rank(W)], the generators of the parabolic subgroup.

intrinsic IHeckeASMod(W::GrpFPCox, I::SeqEnum[RngIntElt]) -> IHkeASMod
{Create the right antispherical module, with basis ^I W (minimal coset representatives for the right
 cosets of W_I in W.}
    I := Sort(Setseq(Seqset(I)));
    require I subset [1..Rank(W)]: "Parabolic subset", I, "should be a subset of", [1..Rank(W)];
    LPolyRing<v> := LaurentSeriesRing(Integers());

    asmod := New(IHkeASMod);
    name := Sprintf("Antispherical module of type %o, parabolic %o", _IHkeCartanName(W), I);
    _FModIHkeInit(~asmod, LPolyRing, name, W);
    asmod`Para := I;
    return asmod;
end intrinsic;

intrinsic 'eq'(ASMod1::IHkeASMod, ASMod2::IHkeASMod) -> BoolElt
{}
    return CoxeterGroup(ASMod1) cmpeq CoxeterGroup(ASMod2) and ASMod1`Para eq ASMod2`Para;
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

intrinsic IHeckeSMod(W::GrpFPCox, I::SeqEnum[RngIntElt]) -> IHkeSMod
{Create the right spherical module, with basis ^I W (minimal coset representatives for the right
 cosets of W_I in W.}
    I := Sort(Setseq(Seqset(I)));
    require I subset [1..Rank(W)]: "Parabolic subset", I, "should be a subset of", [1..Rank(W)];
    LPolyRing<v> := LaurentSeriesRing(Integers());

    smod := New(IHkeSMod);
    name := Sprintf("Spherical module of type %o, parabolic %o", _IHkeCartanName(W), I);
    _FModIHkeInit(~smod, LPolyRing, name, W);
    smod`Para := I;
    return smod;
end intrinsic;

intrinsic 'eq'(SMod1::IHkeSMod, SMod2::IHkeSMod) -> BoolElt
{}
    return CoxeterGroup(SMod1) cmpeq CoxeterGroup(SMod2) and SMod1`Para eq SMod2`Para;
end intrinsic;

intrinsic Parabolic(SMod::IHkeSMod) -> SeqEnum[RngIntElt]
{The parabolic generators.}
    return SMod`Para;
end intrinsic;
