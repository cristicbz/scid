module scid.demos.issue49;

version( demo ):
import scid.demos.common;

/** Issue 49 - Wrong index-slice-assign */
void testIssue49() {
    auto mat = Matrix!double([[1,2,3],[4,5,6],[7,8,9]]);
    auto vec = Vector!double([100, 200]).t;
    mat[2][0..2] = vec;
    enforceMatData( mat, 3, 3, [1.,4,100,2,5,200,3,6,9] );
}
