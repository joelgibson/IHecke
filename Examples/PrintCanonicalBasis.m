// Examples/PrintCanonicalBasis.m
//
// Prints the canonical basis of a given type, in terms of the standard basis. Example usage:
//  $ magma type:=A2 Examples/PrintCanonicalBasis.m
//  $ magma type:=B3 Examples/PrintCanonicalBasis.m
//
// Passing tabular:=true will make it print a table, it may be good to send this table through
// another function to make sure it's properly formatted:
//  $ magma -b type:=A3 tabular:=true Examples/PrintCanonicalBasis.m | column -t -s "\t" -n

AttachSpec("IHecke.spec");
SetColumns(0);
SetQuitOnError(true);

procedure Usage(missingArg)
    print "";
    printf "Error: missing argument %o\n", missingArg;
    print "Usage example: magma type:=B2 Examples/PrintCanonicalBasis.m";
    print "";
    print "Arguments:";
    print "  type    (required) Cartan type, eg A2, B3, F4, ...";
    print "  tabular (optional) Print a tab-separated table, to feed into `column -t`.";
end procedure;

if not assigned type then
    Usage("type");
    quit;
end if;

W := CoxeterGroup(GrpFPCox, type);
HAlg := IHeckeAlgebra(W);
H := IHeckeAlgebraStd(HAlg);
C := IHeckeAlgebraCan(HAlg);

procedure simplePrint()
    printf "Numbering convention for type %o:\n", type;
    CoxeterDiagram(W);
    print "";
    print "Canonical basis:";

    for w in EnumerateCoxeterGroup(W) do
        printf "%o = %o\n", C.w, H ! C.w;
    end for;
end procedure;

procedure tabularPrint()
    // Show a Coxter word in brief form, eg "id" or "121".
    function FmtElt(w)
        return #w eq 0
            select "id"
            else   &cat[IntegerToString(i) : i in Eltseq(w)];
    end function;

    elts := EnumerateCoxeterGroup(W);

    // Header row
    print &*["\t" * FmtElt(w) : w in elts];

    for w in elts do
        printf FmtElt(w);
        Cw := H ! C.w;
        for y in elts do
            coeff := Coefficient(Cw, y);
            if coeff ne 0 then
                printf "\t%o", &*Split(Sprint(coeff), " ");
            else;
                printf "\t";
            end if;
        end for;
        printf "\n";
    end for;
end procedure;

if not assigned tabular then
    simplePrint();
else
    tabularPrint();
end if;

quit;