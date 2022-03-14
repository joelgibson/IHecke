// MuCoeffs.m is an API for reading tables of Mu-coefficients and loading them into the rest of IHecke.
//  1. _LoadMuTable(W) returns (true, table) if a relevant table exists, and false otherwise. The table should be
//     treated as an opaque object.
//  2. _MuCoeffFromTable(table, w) returns an associative array mapping elements x of W to mu(x, w). Zeros do not appear
//     in this map.

// Returns the path to the directory containing this file: returns the directory with no trailing slash.
// Copied from the MAGMA interfaces to Growl, Howl, and Boxcar by Alexander Kasprzyk and Dan Roozemond.
function mucoeff_dir()
    // Magma doesn't provide a way for a package file to discover its
    // path directly -- instead we force an error and parse the error string
    // to extract the path.
    try
        error "Catch me!";
    catch e;
        // First line is: In file "/path/to/file.m", line x, column y. Extract the file path.
        path := Split(e`Position, "\"")[2];

        // Trim the trailing /file.m
        ok, substr := Regexp(".*/", path);
        error if not ok, "Could not parse directory name";
        return substr[1..#substr-1];
    end try;
end function;

MuTableRF := recformat< MuTable : MtrxSprs, WBij : SetIndx>;

function parseMuFile(W, fileContents)
    data := eval fileContents;
    error if data`CoxeterMat ne CoxeterMatrix(W), "Coxeter matrices are incompatible";

    WToIndex := AssociativeArray(W);
    WToIndex[W.0] := 1;
    IndexToW := AssociativeArray(Integers());
    IndexToW[1] := W.0;
    mult := data`ShortLexMultTable;
    for idx in [1..Nrows(mult)] do
        w := IndexToW[idx];

        // Generate any new w elements appearing on this row.
        for s in Support(mult, idx) do
            sw := W.s * w;
            WToIndex[sw] := mult[idx, s];
            IndexToW[mult[idx, s]] := sw;
        end for;
    end for;

    // Write out an indexed set for later reference.
    WBij := {@ IndexToW[i] : i in [1 .. #IndexToW] @};
    return rec<MuTableRF | MuTable:=data`MuTable, WBij:=WBij>;
end function;

intrinsic _MuCoeffFromTable(table::Rec, w::GrpFPCoxElt) -> Assoc
{Retrieve the vector mu(-, w) from the table.}
    mu := AssociativeArray(Parent(w));
    for x in BruhatDescendants(w) do
        mu[x] := 1;
    end for;

    wIdx := Index(table`WBij, w);
    for xIdx in Support(table`MuTable, wIdx) do
        mu[table`WBij[xIdx]] := table`MuTable[wIdx, xIdx];
    end for;

    return mu;
end intrinsic;

intrinsic _LoadMuTable(W::GrpFPCox) -> BoolElt, Assoc
{Load the mu-table for a Coxeter group, returning (true, table) if one exists, and (false, false) otherwise.}
    path := mucoeff_dir() cat "/" cat CartanName(W) cat ".mus";
    ok, io := OpenTest(path, "r");
    if not ok then
        return false, false;
    end if;

    return true, parseMuFile(W, Read(io));
end intrinsic;
