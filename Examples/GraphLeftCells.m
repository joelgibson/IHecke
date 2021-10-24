// Examples/GraphLeftCells.m
//
// Generates a Graphviz file of the left cells. To check the Graphviz output, run
//  $ magma -b type:=G2 prime:=3 Examples/GraphLeftCells.m
//
// To generate the actual picture, run
//  $ magma -b type:=G2 prime:=3 Examples/GraphLeftCells.m | tred | dot -Tpng -o g2p3.png


AttachSpec("IHecke.spec");
SetColumns(0);
SetQuitOnError(true);

procedure Usage(missingArg)
    print "";
    printf "Error: missing argument %o\n", missingArg;
    print "Usage example: magma -b type:=B2 prime:=2 Examples/GraphLeftCells.m";
    print "";
    print "Arguments:";
    print "  type    (required) Cartan type, eg A2, B3, F4, ...";
    print "  prime   (optional) A prime (activates the p-canonical basis) or 0.";
end procedure;

if not assigned type then
    Usage("type");
    quit;
end if;

if not assigned prime then
    prime := 0;
else;
    prime := StringToInteger(prime);
    error if prime ne 0 and not IsPrime(prime),
        "prime must be zero or a prime";
end if;

// Set up the algebra, basis, and cells
W := CoxeterGroup(GrpFPCox, type);
HAlg := IHeckeAlgebra(W);
H := StandardBasis(HAlg);
B := prime eq 0
    select CanonicalBasis(HAlg)
      else PCanonicalBasis(HAlg, type, prime: quiet:=true);
left, right, twoSided := Cells(B);

// Show a Coxter word in brief form, eg "id" or "121".
function FmtElt(w)
    return #w eq 0
        select "id"
        else   &cat[IntegerToString(i) : i in Eltseq(w)];
end function;

// Start outputting Graphviz directions
printf "digraph G {\n";
printf "  graph [ rankdir = \"BT\"; fontsize=11; ];\n";
printf "  node [ fontname=monospace; fontsize=11; ];\n";

// Vertices
comps := Components(left);
for compIdx -> comp in comps do
    printf "  node%o [\n", compIdx;
    printf "    shape = \"record\"\n";
    printf "    label = \"{%o}\"\n",
        Sprintf("LCell #%o", compIdx) * &*["|" * FmtElt(w) : w in Sort(Setseq(comp))];
    printf "  ];\n";
end for;

// Edges
print "";
for edge in GeneratingEdges(left) do
    printf "  node%o -> node%o;\n", Index(comps, edge[1]), Index(comps, edge[2]);
end for;

printf "  label=\"Left cells of %o, p=%o\"\n", type, prime;
printf "  labelloc=\"t\";\n";
printf "}\n";

quit;