module scid.demos.dmatops;

version( demo ):
import scid.demos.common;

/** Test a medley of matrix operations with double elements. */
void dMatOpsTest() {
    alias Matrix!double            dGeMat;
    alias SymmetricMatrix!double   dSyMat;

    auto a = dGeMat( 3, [1.,2,3,4,5,6,7,8,9] );
    auto b = dGeMat( 3, [1.,2,3,4,5,6] );

    dGeMat c = b * a;
    enforceMatData( c, 2, 3, [22,28,49,64,76,100] );
    c[] = c[0 .. 2][ 0 .. 2 ].t * ( (b[][0] - a[1..3][0]).t * eval(c[][0]) ) / 50.;
    enforceMatData( c, 2, 2, [-22,-49,-28,-64] );

// Issue 84
    dSyMat s = c.t*c;

    enforceMatData( s, 2, 2, [2885, 3752, 4880] );

    auto d = eval( s - dSyMat([2800.,3700,4800]) );
    static assert( is( typeof(d) : dSyMat ) );
    enforceMatData( d, 2, 2, [85,52,80] );
    enforce( d[1][0] == 52 );

    auto e = eval( d - b[0..2][1..3]*10 );
    static assert( is( typeof(e) : dGeMat ) );
    enforceMatData( e, 2, 2, [ 55, 12, 2, 20 ] );
}
