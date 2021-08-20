import "Base.m":
    _LaurentPolyRing,
    _IHkeFModInit;

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
{Two Hecke algebras are "equal" if the Coxeter matrices of their Coxeter groups are equal.}
    return algA`CoxeterMatrix eq algB`CoxeterMatrix;
end intrinsic;

intrinsic DefaultBasis(alg::AlgIHke) -> AlgIHkeStd
{}
    return IHeckeAlgebraStd(alg);
end intrinsic;