declare type CelIHke;
declare attributes CelIHke:
    // The basis associated to the cells (some type extending AlgIHkeBase)
    Basis,

    // The CellType is "Left", "Right", or "Two-sided".
    CellType,

    // A component is an indexed set, containing sets of Coxeter group elements.
    // For example, for the two-sided cells of S3 we would have
    // {@ {id}, {s, t, st, ts}, {sts} @}
    Components,

    // An associative array mapping a Coxeter group element to its cell index, eg in the above
    // we would have id => 1, s => 2, t => 2, st => 2, ts => 2, and sts => 3.
    CoxToCompIdx,

    // The quotient digraph, whose transitive closure gives the cell ordering. In the S3 example,
    // the digraph is a path 1 -> 2 -> 3. This quotient is non-canonical, but its transitive
    // closure (or transitive reduction) are canonical.
    QuotGraph;


////////////////////////
// Accessors

intrinsic Print(cel::CelIHke)
{}
    printf "%o %o cells of %o", #cel`Components, cel`CellType, cel`Basis;
end intrinsic;

intrinsic Components(cel::CelIHke) -> SetIndx[SetEnum[GrpFPCox]]
{Return the components (cells) in the decomposition.}
    return cel`Components;
end intrinsic;

intrinsic Component(cel::CelIHke, x::GrpFPCoxElt) -> SetEnum[GrpFPCox]
{Return the cell containing x.}
    xIdx := cel`CoxToCompIdx[x];
    return cel`Components[xIdx];
end intrinsic;

intrinsic CellLe(cel::CelIHke, x::GrpFPCoxElt, y::GrpFPCoxElt) -> BoolElt
{Given a cell decomposition and a pair x, y, determine whether x <= y in the cell order.}
    xIdx := cel`CoxToCompIdx[x];
    yIdx := cel`CoxToCompIdx[y];
    verts := VertexSet(cel`QuotGraph);
    ok := PathExists(verts ! xIdx, verts ! yIdx);
    return ok;
end intrinsic;

// Experimental: returning the edges between cells. Eventually this should return the transitive
// reduction.
intrinsic GeneratingEdges(cel::CelIHke) -> SeqEnum[SetEnum[GrpFPCox]]
{Return a list of edges [comp1, comp2] which generate the cell order.}
    return [ [cel`Components[Index(InitialVertex(edge))], cel`Components[Index(TerminalVertex(edge))] ]
           : edge in Edges(cel`QuotGraph)
           | InitialVertex(edge) ne TerminalVertex(edge) ];
end intrinsic;


/////////////////////
// Creation

// Helper function since we need to iterate over the whole Coxeter group.
intrinsic EnumerateCoxeterGroup(W::GrpFPCox: lengthBound := -1, quiet := false) -> SetIndx[GrpFPCoxElt]
{Enumerate every element of the Coxeter group, returning the results in a list. If lengthBound is
 nonnegative, only elements up to that length will be returned. If lengthBound is negative, all
 elements will be returned. An error will be thrown if W is infinite and lengthBound is negative.
 A warning will be issued if more than 1 million elements are enumerated, which can be suppressed
 with quiet := true.}
    warningLimit := 1000000;
    error if lengthBound lt 0 and not IsCoxeterFinite(CartanName(W)),
        "Cannot enumerate an infinite group.";
    lengthBound := lengthBound lt 0
        select #LongestElement(W)
        else Min(#LongestElement(W), lengthBound);

    elts := {@ W.0 @};
    gens := [W.s : s in [1 .. Rank(W)]];
    frontier := { W.0 };
    maxLength := 0;
    while maxLength lt lengthBound do
        newFrontier := {W|};
        for w in frontier do
            for s in gens do
                ws := w * s;
                if #ws gt #w and ws notin newFrontier then
                    Include(~newFrontier, ws);
                    Include(~elts, ws);
                end if;
            end for;
        end for;
        frontier := newFrontier;
        maxLength +:= 1;
    end while;

    return elts;
end intrinsic;

// Contructor for a single type of cell (left, right, two-sided).
function _CreateCel(Basis, CellType, Welts, digraph)
    // Each cmp is a subset of GrphVert objects rather than Coxeter group elements. I don't know of
    // a better way of converting them back than looking up their index in the vertex set, then
    // using the vertex set again to find the element.
    sccs := {@ comp : comp in StronglyConnectedComponents(digraph) @};

    CoxToCompIdx := AssociativeArray(CoxeterGroup(Basis));
    Components := {@ @};
    for i -> comp in sccs do
        newComp := {};
        for wVert in comp do
            wIdx := Index(Vertices(digraph), wVert);
            w := Welts[wIdx];

            Include(~newComp, w);
            CoxToCompIdx[w] := i;
        end for;
        Include(~Components, newComp);
    end for;

    // I want to use quo<digraph | sccs> here, but that doesn't seem to preserve order, and there
    // doesn't seem to be any other way of getting the quotient map from the quotient graph
    // constructor. We'll (awkwardly) form our own quotient graph.
    edges := {};
    for edge in Edges(digraph) do
        source := Welts[Index(Vertices(digraph), InitialVertex(edge))];
        dest := Welts[Index(Vertices(digraph), TerminalVertex(edge))];
        Include(~edges, [CoxToCompIdx[source], CoxToCompIdx[dest]]);
    end for;
    QuotGraph := Digraph< #Components | edges : SparseRep := true>;

    cel := New(CelIHke);
    cel`Basis := Basis;
    cel`CellType := CellType;
    cel`Components := Components;
    cel`CoxToCompIdx := CoxToCompIdx;
    cel`QuotGraph := QuotGraph;
    return cel;
end function;


// User-callable, constructs all cells.
intrinsic Cells(B::AlgIHkeBase) -> Digraph, Digraph, Digraph
{Return the left, right, and two-sided cells of B. Each collection of cells is returned as a
 directed acyclic graph, whose vertices are sets of Coxeter group elements, forming a partition of
 the Coxeter group. The transitive closure of the graph gives the cell order.}
    W := CoxeterGroup(B);
    Welts := EnumerateCoxeterGroup(W);
    leftEdges := {};
    rightEdges := {};
    for w in Welts do
        for s in [1..Rank(W)] do
            left := B ! (B.s * B.w);
            for u in Support(left) do
                Include(~leftEdges, [w, u]);
            end for;

            right := B ! (B.w * B.s);
            for u in Support(right) do
                Include(~rightEdges, [w, u]);
            end for;
        end for;
    end for;

    edgeSets := [leftEdges, rightEdges, leftEdges join rightEdges];
    digraphs := [Digraph<Welts | edges : SparseRep := true> : edges in edgeSets];
    return
        _CreateCel(B, "Left", Welts, digraphs[1]),
        _CreateCel(B, "Right", Welts, digraphs[2]),
        _CreateCel(B, "Two-sided", Welts, digraphs[3]);
end intrinsic;