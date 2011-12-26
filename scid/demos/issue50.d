module scid.demos.issue50;

version( demo ):
import scid.demos.common;

/** Issue 50 - Matrix slice-slice-assign ends up transposed	*/
void testIssue50() {
    auto a = Matrix!double([[1.0, 2], [4.0, 5]]);
    auto b = Matrix!double(3, 3);
    b[0..2][0..2] = a;
    enforce( b[0,0] == 1. && b[0,1] == 2. && b[1,0] == 4. && b[1,1] == 5. );
}
