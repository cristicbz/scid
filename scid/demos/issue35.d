module scid.demos.issue35;

version( demo ):
import scid.demos.common;

/** Issue 35 - Assignment should return value for chaining */
void testIssue35() {
    auto mat = Matrix!double(5, 5);
    auto vec = Vector!double([1.,2,3,4]);
    mat[2, 2] = mat[3, 3] = 5;
    enforce( mat[2, 2] == 5 && mat[3, 3] == 5 );
    mat[1, 1] = mat[2, 2] += 3;
    enforce( mat[1, 1] == 8 && mat[2, 2] == 8);

    vec[ 0 ] = vec[ 1 ] = 10.;
    enforce( vec[0] == 10. && vec[1] == 10. );
    vec[ 2 ] = vec[ 3 ] += 5.0;
    enforce( vec[2] == 9.0 && vec[ 3 ] == 9.0 );
}
