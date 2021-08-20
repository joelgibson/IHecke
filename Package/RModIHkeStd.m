import "Base.m":
    _LaurentPolyRing,
    _v;
import "EltIHke.m":
    _EltIHkeConstruct,
    _AddScaled,
    _RemoveZeros,
    _AddScaledTerm;
import "RModIHke.m":
    _RModIHeckeBaseInit;

declare type RModIHkeStd[EltIHke]: RModIHkeBase;

intrinsic IHeckeRModuleStd(alg::AlgIHke, para::SeqEnum[RngIntElt]) -> RModIHkeStd
{Standard basis of the antispherical right module.}
    // TODO: Cache the basis on the algebra object.
    basis := New(RModIHkeStd);
    _RModIHeckeBaseInit(~basis, alg, "aH", "Standard basis", para);
    return basis;
end intrinsic;

intrinsic 'eq'(H1::RModIHkeStd, H2::RModIHkeStd) -> BoolElt
{}
    return Parent(H1) eq Parent(H2) and Parabolic(H1) eq Parabolic(H2);
end intrinsic;


//////////////////////////////////
// Right action

// Right-multiply an element in the standard basis by H(s). Precondition: elt is only supported
// over I-minimal elements.
function _RightMultStdGen(elt, s)
    assert ISA(Type(Parent(elt)), RModIHkeStd);

    W := CoxeterGroup(Parent(elt));
    I := Parabolic(Parent(elt));
    terms := AssociativeArray();
    for w -> coeff in elt`Terms do
        ws := w * (W.s);
        if not IsMinimal(I, ws) then
            _AddScaledTerm(~terms, w, -_v*coeff);
        elif #w lt #ws then
            _AddScaledTerm(~terms, ws, coeff);
        else
            _AddScaledTerm(~terms, ws, coeff);
            _AddScaledTerm(~terms, w, (_v^-1 - _v) * coeff);
        end if;
    end for;
    _RemoveZeros(~terms);
    return _EltIHkeConstruct(Parent(elt), terms);
end function;

intrinsic _IHkeProtAction(M::RModIHkeStd, eltM::EltIHke, H::AlgIHkeStd, eltH::EltIHke) -> EltIHke
{Right action of the standard basis of the Hecke algebra on the standard basis of the module.}
    require Parent(M) eq Parent(H): "Algebra objects must be equal";

    acc := AssociativeArray();
    for w -> coeff in eltH`Terms do
        piece := eltM;
        for s in Eltseq(w) do
            piece := _RightMultStdGen(piece, s);
        end for;
        _AddScaled(~acc, piece`Terms, coeff);
    end for;
    _RemoveZeros(~acc);
    return _EltIHkeConstruct(M, acc);
end intrinsic;