module scid.demos.issue61;

version( demo ):
import scid.demos.common;

/** Issue 61 - Transposing eagerly copies */
void testIssue61() {
    auto mat1 = Matrix!double([ [1.,2.], [3., 4.] ]);
    auto mat2 = eval( mat1.t );
    Matrix!double mat3;
    mat3[] = mat2.t;

    enforce( mat1.cdata == mat2.cdata );
    enforce( mat2.cdata == mat3.cdata );
}
