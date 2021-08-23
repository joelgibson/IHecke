// This example shows how to define a new basis of the Hecke algebra, the Kazhdan-Lusztig-normalised
// basis. Its relationship to our usual standard basis is that H(x) = v^#x T(x).
//
// Since this file defines intrinsics, it needs to be loaded through a Spec file (KLStandard.spec).

// First a new type `AlgIHkeKLStd` needs to be created, inheriting from the abstract type BasisIHke.
// Inheriting from BasisIHke means that many operations are already defined, such as printing,
// accessing the base ring, the equality predicate, constructing elements, and so on (see
// Package/BasisIHke.m for the full list). We have also declared an element type of `EltIHke`, which
// tells Magma that "elements" of this structure are of type `EltIHke`.
declare type AlgIHkeKLStd[EltIHke]: BasisIHke;

// Next, define a function to create an instance of the basis from an instance of a Hecke algebra.
// This is mostly just creating a type - however we cache the basis object so we don't have to
// create it again later (this is not so important here, but important if bases are carrying
// expensive-to-compute cached information that should not be recomputed).
intrinsic IHeckeAlgebraKLStd(HAlg::AlgIHke) -> AlgIHkeKLStd
{The Kazhdan-Lusztig-normalised standard basis of the Hecke algebra.}
    // If a basis object has already been cached, return it.
    if IsDefined(HAlg`BasisCache, AlgIHkeKLStd) then
        return HAlg`BasisCache[AlgIHkeKLStd];
    end if;

    // Otherwise, create a new basis object. It is passed to _BasisIHkeInit(...), which initialises
    // some basic parameters.
    T := New(AlgIHkeKLStd);
    _BasisIHkeInit(~T, HAlg, "T", "KL-standard basis");

    // Save the new basis in the cache, and return it.
    HAlg`BasisCache[AlgIHkeKLStd] := T;
    return T;
end intrinsic;

// At this point, expressions like T.0 + T.1 - T.[2, 3] and so on will work just fine. To be able
// to change into other bases, we need to define two basis conversion functions: one from this new
// basis into the standard basis, and one from the standard basis into this new basis.
//
// We will define the conversion from T to the standard basis H first. To do this we only need to
// define an overload for the intrinsic `_IHkeProtToBasis` for the correct types: the '*' operation
// will automatically call this.
intrinsic _IHkeProtToBasis(H::AlgIHkeStd, T::AlgIHkeKLStd, w::GrpFPCoxElt) -> EltIHke
{Convert from the KL-standard basis to the standard basis.}
    // To get a hold of the base ring, either H or T can be used (when this function is invoked
    // automatically, H and T will be coming from the same parent Hecke algebra, and hence hold
    // equal base rings and Coxeter groups). The '.1' after the base ring retrieves the first (and
    // only) generator v.
    v := BaseRing(H).1;

    // The length of a Coxeter group element in Magma is #w.
    return H.w * v^(-#w);
end intrinsic;

// The basis conversion in the other direction is just as simple.
intrinsic _IHkeProtToBasis(T::AlgIHkeKLStd, H::AlgIHkeStd, w::GrpFPCoxElt) -> EltIHke
{Convert from the standard basis to the KL-standard basis.}
    v := BaseRing(H).1;
    return T.w * v^(#w);
end intrinsic;

// And that's all there is to it. Operations like multiplication and the bar involution will now
// automatically work on this basis (they will transform into the standard basis, perform the
// operation there, and transform back). If these automatic operations are too inefficient, each
// operation can be overloaded for each basis or pair of bases à la carte.