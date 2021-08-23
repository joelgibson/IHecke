import "EltIHke.m": _AddScaledTerm, _RemoveZeros;

///////////////////////////////////////
// The group algebra of a Coxeter group

declare type IHkeGrpAlg: FModIHke;

intrinsic IHeckeGroupAlgebra(W::GrpFPCox) -> IHkeGrpAlg
{Create a new group algebra object.}
    GAlg := New(IHkeGrpAlg);
    name := Sprintf("Group algebra of type %o", CartanName(W));
    _FModIHkeInit(~GAlg, Integers(), name, W);
    return GAlg;
end intrinsic;

intrinsic 'eq'(GAlg1::IHkeGrpAlg, GAlg2::IHkeGrpAlg) -> BoolElt
{}
    return GAlg1`CoxeterMatrix eq GAlg2`CoxeterMatrix;
end intrinsic;

intrinsic DefaultBasis(GAlg::IHkeGrpAlg) -> IHkeGrpAlgStd
{}
    return StandardBasis(GAlg);
end intrinsic;


//////////////////////////////////////
// Standard basis of the group algebra

declare type IHkeGrpAlgStd[EltIHke]: BasisIHke;

intrinsic StandardBasis(GAlg::IHkeGrpAlg) -> IHkeGrpAlgStd
{}
    G := New(IHkeGrpAlgStd);
    _BasisIHkeInit(~G, GAlg, "G", "Standard basis");
    return G;
end intrinsic;

intrinsic _IHkeProtUnit(G::IHkeGrpAlgStd) -> EltIHke
{}
    return G.0;
end intrinsic;

intrinsic _IHkeProtMult(G1::IHkeGrpAlgStd, elt1::EltIHke, G2::IHkeGrpAlgStd, elt2::EltIHke) -> EltIHke
{Multiplication in the standard basis of the group algebra.}
    terms := AssociativeArray(CoxeterGroup(G1));
    for w1 -> coeff1 in elt1`Terms do
        for w2 -> coeff2 in elt2`Terms do
            _AddScaledTerm(~terms, w1 * w2, coeff1 * coeff2);
        end for;
    end for;
    _RemoveZeros(~terms);
    return EltIHkeConstruct(G1, terms);
end intrinsic;

intrinsic _IHkeNaturalMap(G::IHkeGrpAlgStd, H::IHkeAlgStd, eltH::EltIHke) -> EltIHke
{The natural map from the Hecke algebra to the group algebra by specialisation at v=1.}
    terms := AssociativeArray(CoxeterGroup(G));
    for w -> coeff in eltH`Terms do
        _AddScaledTerm(~terms, w, Evaluate(coeff, 1));
    end for;
    return EltIHkeConstruct(G, terms);
end intrinsic;

intrinsic _IHkeProtMult(G::IHkeGrpAlgStd, eltG::EltIHke, H::IHkeAlgStd, eltH::EltIHke) -> EltIHke
{The right action of the Hecke algebra on the group algebra by specialisation at v=1.}
    return eltG * _IHkeNaturalMap(G, H, eltH);
end intrinsic;

intrinsic _IHkeProtMult(H::IHkeAlgStd, eltH::EltIHke, G::IHkeGrpAlgStd, eltG::EltIHke) -> EltIHke
{The left action of the Hecke algebra on the group algebra by specialisation at v=1.}
    return _IHkeNaturalMap(G, H, eltH) * eltG;
end intrinsic;