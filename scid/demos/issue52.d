module scid.demos.issue52;

version( demo ):
import scid.demos.common;

/** Issue 52 - Inverting Scalars */
void testIssue52() {
    // Scalar literals
    float   f = 2.0f, invf = inv( f );
    double  d = 4.0,  invd = inv( d );
    cdouble z = 2.0 + 2.0i, invz = inv( z );
    enforce( invf == 0.5f );
    enforce( invd == 0.25 );
    enforce( invz == 0.25 - 0.25i );

    // Scalar expression
    double invdot = eval( inv( externalVectorView!(VectorType.Row)( [1.,2.] ) * externalVectorView( [4.,8.] ) ) );
    enforce( invdot == 0.05 );
}
