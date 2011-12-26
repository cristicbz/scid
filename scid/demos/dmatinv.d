module scid.demos.dmatinv;

version( demo ):
import scid.demos.common;

/** Test inversions on matrices of doubles. */
void dMatInvTest() {
    auto x = Matrix!double([ [ 1, 2, 3 ], [  1, 1, 1 ], [ 2, -1, 1 ] ]);
    auto y = Matrix!double([ [ 1, 2, 0 ], [ -2, 3, 4 ], [ 0,  2, 1 ] ]);

    // Thank you Octave...
    enforceMatData( eval(inv(x)*y), 3, 3,     [  -2.4,  -2.2,   2.6,   2.6,   1.8,  -1.4,   4.2,   3.6,  -3.8 ] );
//    enforceMatData( eval(y*inv(x)), 3, 3,     [  -0.8,   2.6,   0.2,   3.0,  -3.0,   1.0,  -0.6,  -0.8,  -0.6 ] );  // Issue 82
//    enforceMatData( eval(inv(x.t)*y), 3, 3,   [   0.0,  -1.0,   1.0,  -0.2,   3.0,  -0.4,  -0.2,   3.0,  -1.4 ] );  // Issue 83
    enforceMatData( eval(y*inv(x.t)), 3, 3,   [   1.6,   4.6,   2.2,   1.8,   1.8,   1.6,  -1.4,  -3.4,  -1.8 ] );
    enforceMatData( eval(inv(y)*x), 3, 3,     [  -9.0,   5.0,  -8.0,  20.0,  -9.0,  17.0,   9.0,  -3.0,   7.0 ] );
//    enforceMatData( eval(x*inv(y)), 3, 3,     [  13.0,   7.0,  16.0,   6.0,   3.0,   7.0, -21.0, -11.0, -27.0 ] );  // Issue 82
    enforceMatData( eval(inv(y.t)*x), 3, 3,   [  11.0,   5.0, -18.0,   4.0,   1.0,  -5.0,  17.0,   7.0, -27.0 ] );
    enforceMatData( eval(x*inv(y.t)), 3, 3,   [ -15.0,  -1.0,   0.0,   8.0,   1.0,   1.0, -13.0,  -1.0,  -1.0 ] );

    enforceMatData( eval(x.t*inv(y.t)), 3, 3, [  -9.0,  20.0,   9.0,   5.0,  -9.0,  -3.0,  -8.0,  17.0,   7.0 ] );
    enforceMatData( eval(y.t*inv(x.t)), 3, 3, [  -2.4,   2.6,   4.2,  -2.2,   1.8,   3.6,   2.6,  -1.4,  -3.8 ] );  // Issue 82
}
