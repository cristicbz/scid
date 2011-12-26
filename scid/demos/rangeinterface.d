module scid.demos.rangeinterface;

version( demo ):
import scid.demos.common;

/** Using vectors and matrices and ranges. */
void rangeInterface() {
    writeln();
    writeln( "======================== Range Interface ========================" );

    // InputRange. Using foreach with vectors
    auto v = Vector!double([1.0, 2.0, 3.0]);
    auto sum = 0.;
    foreach( e ; v ) sum += e;
    writeln( "The sum of ", v.toString(), "'s elements is ", sum );

    // Matrices are iterated by major subvector (e.g. columns for column-major matrices). Importantly, the elements
    // in the iterations are views. Therefore changing an element affects the matrix, no matter if ref is used or not.
    auto rowMat = Matrix!(double, StorageOrder.RowMajor)([ [1.,2.,3.], [4.,5.,6.], [7.,8.,9.] ]);
    auto colMat = Matrix!double( [ [1.,2.,3.], [4.,5.,6.], [7.,8.,9.] ] );

    uint i = 0;
    writeln( "Row major matrix: " );
    foreach( r ; rowMat ) {
        writeln( "Row ", i, ": ", r.toString() );

        if( i == 0 )      enforce( r == [1.,2.,3.] );
        else if( i == 1 ) enforce( r == [4.,5.,6.] );
        else if( i == 2 ) enforce( r == [7.,8.,9.] );

        i ++;
    }
    writeln();

    i = 0;
    writeln( "Column major matrix: " );
    foreach( c ; colMat ) {
        writeln( "Column ", i, ": ", c.toString() );

        if( i == 0 )      enforce( c == [1.,4.,7.] );
        else if( i == 1 ) enforce( c == [2.,5.,8.] );
        else if( i == 2 ) enforce( c == [3.,6.,9.] );

        i ++;
    }
    writeln();
}
