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

function parseMuFile(W, fileContents)
    lines := Split(fileContents, "\n");

    // The first line gives all the reduced expressions, in index order, with indices starting from 0.
    index := [W ! word : word in eval lines[1]];

    // The following lines are in index order, and map elements to arrays of pairs, so that on the
    // line corresponding to w for example, we have
    // [x, mu(x, w), y, mu(y, w), ...] for each nonzero value of mu.
    mus := AssociativeArray(W);
    for i -> pairs in [eval line : line in lines[2..#lines]] do
        w := index[i];
        mus[w] := AssociativeArray(W);
        for j in [1 .. #pairs by 2] do
            x := index[pairs[j] + 1];
            mus[w][x] := pairs[j+1];
        end for;
    end for;

    return mus;
end function;

intrinsic _LoadMuCoefficients(W::GrpFPCox) -> BoolElt, Assoc
{Attempt to load the mu-coefficients for a given Coxeter group from the database in IHecke, returning
 a status code (false: not found, true: available), and an associative array mapping group elements
 to more associative arrays, so that mu(x, y) is at [y][x] in the resulting data structure.}
    // Note that the Cartan name of type Cn is Bn.
    path := mucoeff_dir() cat "/" cat CartanName(W) cat ".mus";
    ok, io := OpenTest(path, "r");
    if not ok then
        return false, AssociativeArray(W);
    end if;

    return true, parseMuFile(W, Read(io));
end intrinsic;

