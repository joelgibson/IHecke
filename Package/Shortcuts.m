// This file defines some "shortcut methods" for defining an algebra or module, along with important
// bases, all at once (to simplify interpreter sessions).

intrinsic ShortcutIHeckeAlgebra(W::GrpFPCox) -> IHkeAlg, GrpFPCox, IHkeAlgStd, IHkeAlgCan
{Shortcut for constructing the Hecke algebra, group, standard, and canonical bases, from a group.}
    HAlg := IHeckeAlgebra(W);
    return HAlg, W, IHeckeAlgebraStd(HAlg), IHeckeAlgebraCan(HAlg);
end intrinsic;

intrinsic ShortcutIHeckeAlgebra(type::MonStgElt) -> IHkeAlg, GrpFPCox, IHkeAlgStd, IHkeAlgCan
{Shortcut for constructing the Hecke algebra, group, standard, and canonical bases, from a type.}
    return ShortcutIHeckeAlgebra(CoxeterGroup(GrpFPCox, type));
end intrinsic;

intrinsic ShortcutIHeckeAntiSpherical(HAlg::IHkeAlg, I::SeqEnum[RngIntElt]) -> IHkeASMod, ASModIHkeStd, ASModIHkeCan
{Shortcut for constructing the antispherical module, standard, and canonical bases from a Hecke algebra.}
    ASMod := IHeckeAntiSpherical(HAlg, I);
    return ASMod, IHeckeAntiSphericalStd(ASMod), IHeckeAntiSphericalCan(ASMod);
end intrinsic;

intrinsic ShortcutIHeckeSpherical(HAlg::IHkeAlg, I::SeqEnum[RngIntElt]) -> IHkeSMod, SModIHkeStd, SModIHkeCan
{Shortcut for constructing the spherical module, standard, and canonical bases from a Hecke algebra.}
    SMod := IHeckeSpherical(HAlg, I);
    return SMod, IHeckeSphericalStd(SMod), IHeckeSphericalCan(SMod);
end intrinsic;