module scid.demos.datainterface;

version( demo ):
import scid.demos.common;

/** Directly accessing the data of vectors and matrices. */
void dataInterface() {
    writeln();
    writeln( "========================= Data Interface =========================" );

    // Most storages provide data & cdata methods which allow access to the raw memory that can be passed to BLAS
    // or used by custom function.
    // cdata() - returns a const pointer to a memory block which, due to copy-on-write, might be shared.
    // data()  - returns a mutable pointer to the memory block. It first ensures that the memory is not shared.

    // By printing the memory of the two following matrices we can see the difference between the storage
    // orders:
    auto rowMat = Matrix!(double, StorageOrder.RowMajor)([ [1.,2.,3.], [4.,5.,6.], [7.,8.,9.] ]);
    auto colMat = Matrix!double( [ [1.,2.,3.], [4.,5.,6.], [7.,8.,9.] ] );

    writeln( "Row major data   : ", rowMat.cdata[ 0 .. 9 ] );
    writeln( "Column major data: ", colMat.cdata[ 0 .. 9 ] );

    // Assigning colMat to a new matrix will cause the two matrices to share the data
    auto otherMat = colMat;
    enforce( otherMat.cdata == colMat.cdata );

    // Calling data will cause the memory to be copied though:
    enforce( otherMat.data  != colMat.data );
    enforce( otherMat.cdata != colMat.cdata );
}
