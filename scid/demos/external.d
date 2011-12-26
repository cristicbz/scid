module scid.demos.external;

version( demo ):
import scid.demos.common;
import scid.internal.regionallocator;

/** Creating views of existing data that allow it to be treated as a matrix or vector. */
void externalViews() {
    auto alloc = newRegionAllocator();
    double array[] = [1.,4,4,5];

    // ===== ExternalMatrixView =====
    // Provides a matrix view to an external (not managed by scid) memory area.

    // - C-tor taking the major dimension (no. of columns for ColumnMajor matrices) and an array which it uses as memory
    auto arrayMat = ExternalMatrixView!double( 2, array );
    arrayMat[ 0, 0 ] = 2.0;
    enforce( array[ 0 ] == 2.0 );

    // - C-tor taking an initializer and an allocator. Copies the initializer into an array allocated with alloc.
    auto initMat = ExternalMatrixView!double( [[1.0, 4], [4.0, 5]], alloc );

    enforce( initMat[ 0, 0 ] == 1.0 && initMat[ 1, 1 ] == 5.0 );

    // - C-tor taking matrix dimensions and an allocator. Allocates new memory using alloc.
    auto sizeMat = ExternalMatrixView!double( 2, 2, alloc );

    enforce( sizeMat.rows == 2 && sizeMat.columns == 2 );


    // ===== ExternalVectorView =====
    // Provides a vector view to an external (not managed by scid) memory area.
    // - Use array for memory.
    auto arrayVec = ExternalVectorView!double( array );

    // - Create a copy of an array using an allocator
    auto initVec = ExternalVectorView!double( [1.,2,3,4], alloc );

    // - Allocate a new array of given size using an allocator
    auto sizeVec = ExternalVectorView!double( 4, alloc );

    // NOTE: As their names suggest, these two data types are views and, therefore, subject to the aliasing problem
    //       which is not yet solved. This means that, to be on the safe side, whenever using a view on the left
    //       hand side of an assignemnt, don't do:
    //
    //          mat[] += mat*inv(mat);
    //
    //       Do:
    //
    //          mat[] += eval( mat*inv(mat) [, alloc] );
    //
    //       This'll allocate a temporary first and then copy it to the left-hand side.

    // ===== Saving results into arrays =====
    auto mat = Matrix!double([ [1.,2], [3.,4] ]);
    auto vec = Vector!double([ 1., 2. ] );

    // Evaluate the expression and save the result into array.
    auto v = eval(mat*vec, array[ 0 .. 2 ]);

    // The array is now equal to the result of the operation.
    enforce( array[ 0 .. 2 ] == [ 5., 11. ] );

    // As an added bonus, eval() returns an external view to the array.
    static assert( is( typeof( v ) : ExternalVectorView!double ) );
    v[] *= 2.0;
    enforce( array[ 0 .. 2 ] == [ 10., 22. ] );
}

