// This file defines some "shortcut methods" for defining an algebra or module, along with important
// bases, all at once (to simplify interpreter sessions).

intrinsic ShortcutIHeckeAlgebra(W::GrpFPCox) -> IHkeAlg, IHkeAlgStd, IHkeAlgCan
{Shortcut for constructing the Hecke algebra, group, standard, and canonical bases from a group.}
    HAlg := IHeckeAlgebra(W);
    return HAlg, StandardBasis(HAlg), CanonicalBasis(HAlg);
end intrinsic;

intrinsic ShortcutIHeckeASMod(W::GrpFPCox, I::SeqEnum[RngIntElt]) -> IHkeASMod, ASModIHkeStd, ASModIHkeCan
{Shortcut for constructing the antispherical module, standard, and canonical bases from a group.}
    ASMod := IHeckeASMod(W, I);
    return ASMod, StandardBasis(ASMod), CanonicalBasis(ASMod);
end intrinsic;

intrinsic ShortcutIHeckeSMod(W::GrpFPCox, I::SeqEnum[RngIntElt]) -> IHkeSMod, SModIHkeStd, SModIHkeCan
{Shortcut for constructing the spherical module, standard, and canonical bases from a group.}
    SMod := IHeckeSMod(W, I);
    return SMod, StandardBasis(SMod), CanonicalBasis(SMod);
end intrinsic;