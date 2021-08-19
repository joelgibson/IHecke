// The base ring we work over for the Hecke algebra is the Laurent polynomials.
// In Magma, a built-in structure that does this is the Laurent series.
_LaurentPolyRing<v> := LaurentSeriesRing(Integers());
AssignNames(~_LaurentPolyRing, ["v"]);
_v := _LaurentPolyRing.1;

intrinsic IHeckeVersion() -> MonStgElt
{Report version information for IHecke.}
    return "IHecke version 2021-08-18";
end intrinsic;