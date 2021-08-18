import "EltIHke.m": _LaurentPolyRing;

intrinsic IHeckeVersion() -> MonStgElt
{Report version information for IHecke.}
    return "IHecke version 2021-08-18";
end intrinsic;

declare type AlgIHke;
declare attributes AlgIHke:
    // A Coxeter group of finitely-presented type, i.e. a GrpFPCox.
    CoxeterGroup,

    // The Coxeter matrix of CoxeterGroup, used for equality calculations.
    CoxeterMatrix,

    // An associative array with keys type names like AlgIHkeStd, AlgIHkeCan, and so on.
    BasisCache;


///////////////////////
// Accessors

intrinsic Print(alg::AlgIHke)
{A short human-readable description of the algebra.}
    printf "Iwahori-Hecke algebra of type %o", CartanName(alg`CoxeterGroup);
end intrinsic;

intrinsic CoxeterGroup(alg::AlgIHke) -> GrpFPCox
{The underlying Coxeter group (actually Coxeter system, since a GrpFPCox comes with generators).}
    return alg`CoxeterGroup;
end intrinsic;

intrinsic BaseRing(alg::AlgIHke) -> Rng
{The Laurent series ring (we use only the polynomial part).}
    return _LaurentPolyRing;
end intrinsic;

intrinsic 'eq'(algA::AlgIHke, algB::AlgIHke) -> BoolElt
{Two Hecke algebras are "equal" if the Coxeter matrices of their Coxeter groups are equal.}
    return algA`CoxeterMatrix eq algB`CoxeterMatrix;
end intrinsic;


///////////////////////
// Creation

intrinsic IHeckeAlgebra(W::GrpFPCox) -> AlgIHke
{Create a new Hecke algebra.}
    alg := New(AlgIHke);
    alg`CoxeterGroup := W;
    alg`CoxeterMatrix := CoxeterMatrix(W);
    alg`BasisCache := AssociativeArray();
    return alg;
end intrinsic;