// This file defines some "shortcut methods" for defining an algebra or module, along with important
// bases, all at once (to simplify interpreter sessions).

intrinsic ShortcutIHeckeAlgebra(W::GrpFPCox) -> AlgIHke, GrpFPCox, AlgIHkeStd, AlgIHkeCan
{Shortcut for constructing the Hecke algebra, group, standard, and canonical bases, from a group.}
    HAlg := IHeckeAlgebra(W);
    return HAlg, W, IHeckeAlgebraStd(HAlg), IHeckeAlgebraCan(HAlg);
end intrinsic;

intrinsic ShortcutIHeckeAlgebra(type::MonStgElt) -> AlgIHke, GrpFPCox, AlgIHkeStd, AlgIHkeCan
{Shortcut for constructing the Hecke algebra, group, standard, and canonical bases, from a type.}
    return ShortcutIHeckeAlgebra(CoxeterGroup(GrpFPCox, type));
end intrinsic;

// TODO: Rename ASphIHke, ASModIHkeStd, etc, to be consistent.
intrinsic ShortcutIHeckeAntiSpherical(HAlg::AlgIHke, I::SeqEnum[RngIntElt]) -> ASphIHke, ASModIHkeStd, ASModIHkeCan
{Shortcut for constructing the antispherical module, standard, and canonical bases from a Hecke algebra.}
    ASMod := IHeckeAntiSpherical(HAlg, I);
    return ASMod, IHeckeAntiSphericalStd(ASMod), IHeckeAntiSphericalCan(ASMod);
end intrinsic;

intrinsic ShortcutIHeckeSpherical(HAlg::AlgIHke, I::SeqEnum[RngIntElt]) -> SphIHke, SModIHkeStd, SModIHkeCan
{Shortcut for constructing the spherical module, standard, and canonical bases from a Hecke algebra.}
    SMod := IHeckeSpherical(HAlg, I);
    return SMod, IHeckeSphericalStd(SMod), IHeckeSphericalCan(SMod);
end intrinsic;