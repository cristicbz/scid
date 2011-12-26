module scid.demos.issue51;

version( demo ):
import scid.demos.common;

/** Issue 51 - Invalid multiplication between ColumnVector and Matrix */
void testIssue51() {
    // Testing row vectors, too.
    auto col = vector( [1.0, 2, 3] );
    auto row = vector( [4.0, 5, 6] );
    auto mat = matrix( [[4.0, 5, 6]] );

    auto rowRes = eval( col * row.t );
    auto matRes = eval( col * mat );

    foreach( res; [rowRes, matRes] ) {
        assert(	res[0, 0] == 4 );
        assert(	res[0, 1] == 5 );
        assert(	res[0, 2] == 6 );
        assert(	res[1, 0] == 8 );
        assert(	res[1, 1] == 10 );
        assert(	res[1, 2] == 12 );
        assert(	res[2, 0] == 12 );
        assert(	res[2, 1] == 15 );
        assert(	res[2, 2] == 18 );
    }
}
