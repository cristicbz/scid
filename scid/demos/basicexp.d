module scid.demos.basicexp;

version( demo ):
import scid.demos.common;

/** Syntax for basic expressions. */
void basicExpressions() {
    writeln();
    writeln( "======================= Basic Expressions =======================" );
    // The simplest constructors take a built-in array.
    auto mat = Matrix!double( [ [1., 2., 3.], [2., 3., 4.] ] );
    auto vec = Vector!double( [1.,2.,3.,4.,5.,6.] );

    // Printing matrices and vectors
    writeln( "vec = ", vec.toString() ); // Using toString
    writeln( "mat = " );
    writeln( mat.pretty );               // pretty prints rows on multiple lines

    writeln();

    // Slicing vectors and matrices
    auto vecSlice = vec[ 1 .. 3 ];
    auto matSlice = mat[ 0 .. 2 ][ 1 .. 3 ];

    // RowVector times Matrix. The t property transposes vectors and matrices.
    auto w = eval( vecSlice.t * mat );
    writeln( "Expr 1: ", vecSlice.toString(), " * ", mat.toString(), " = ", w.toString() );

    enforce( w == [8.0, 13.0, 18.0] );

    // More complicated expression. mat[0][0..2] gets a view of a slice of the first row. The .t is neccessary since
    // vec is a column vector while the matrix slice is a row vector.
    vecSlice[] = vec[ 0 .. 2 ] * 5.0 - mat[0][0..2].t;
    writeln( "Expr 2: ", vec[0 .. 2].toString(), " * 5 - ", mat[0][0..2].toString(), " = ", vecSlice.toString() );

    enforce( vecSlice == [4.0, 8.0] );

    // (functionality disabled for the moment)
    // One can use array literals as vector literals most of the time:
    // double x = eval( [2.0, -1.0].t * vecSlice );
    // writeln( "Expr 3: [2.0, -1.0].t * ", vecSlice.toString, " = ", x );
    // enforce( x == 0.0 );
}
