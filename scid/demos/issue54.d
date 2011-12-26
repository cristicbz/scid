module scid.demos.issue54;

version( demo ):
import scid.demos.common;

/** Isuue 54 - Matrix slices have copy semantics, while row/column slices have view semantics */
void testIssue54() {
    auto mat   = Matrix!double([ [1.,2.], [3., 4.] ]);
    auto row   = mat[ 0 ][ 0 .. 2 ];
    auto col   = mat[ 0 .. 2 ][ 0 ];
    auto slice = mat[ 0 .. 1 ][ 0 .. 1 ];

    // Copy semantics of matrix slices:
    slice[ 0, 0 ] = 10.0;
    enforce( mat[ 0, 0 ] = 1. );

    // View semantics of row and column slices:
    row[ 0 ] = 100.;
    enforce( mat[ 0, 0 ] == 1. && col[ 0 ] == 1. );

    col[ 0 ] = 200.;
    enforce( mat[ 0, 0 ] == 1. && row[ 0 ] == 100. );
}
