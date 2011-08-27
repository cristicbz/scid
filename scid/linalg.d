/**
This module contains linear algebra functions beyond basic expressions,
such as matrix factorizations.  This module also contains destructive
versions of some functions.  These are named someFunctionDestructive,
and use their input arguments as scratch space.
*/
module scid.linalg;

import std.string;

import scid.internal.regionallocator;
import scid.lapack;
import scid.matrix;
import scid.common.traits;

/**
Performs Cholesky decomposition of a Hermitian, positive definite matrix.  
Assuming the input matrix is named A, Cholesky decomposition has the 
following forms:

mat = L * L.t if L is a lower triangular matrix.
mat = U.t * U if U is an upper triangular matrix.

mat must be symmetric, but it may be stored redundantly using general matrix
storage.  If general matrix storage is used, symmetry is not checked in
the interest of efficiency.

Returns:  

The matrix U that would be multiplied by its conjugate transpose as
shown above to form mat.
*/
auto cholesky( M )( M mat )
if( isMatrix!M ) {  
    alias M.ElementType E;
    auto ret = TriangularMatrix!E( mat.rows, null );
    choleskyImpl( mat, ret );
    return ret;
}

/// Ditto
auto cholesky( M, Allocator )( M mat, Allocator alloc ) 
if( isMatrix!M && isAllocator!Allocator ) {
    alias M.ElementType E;
    alias TriangularMatrix!(E).Temporary Temp;
    auto ret = Temp( mat.rows, alloc );
    choleskyImpl( mat, ret );
    return ret;
}

     
private void choleskyImpl( M1, M2 )( ref M1 mat, ref M2 ret ) {
    alias M1.ElementType E;
    
    auto alloc = newRegionAllocator();
    auto temp = ExternalMatrixView!E( mat.rows, mat.rows, alloc );
    
    // Copy only the upper triangle of mat
    foreach( col; 0..mat.columns ) {
        foreach( row; 0..col + 1 ) {
            temp[ row, col ] = mat[ row, col ];
        }
    }        
    
    int info;    
    lapack.potrf!( 'U' )
        ( temp.length, temp.data, temp.rows, info );
    
    // TODO:  Establish a consistent standard for error handling.
    // For now, just return a matrix of NaNs.
    if( info > 0 ) {
        ret.data[ 0..ret.length ] = E.init;
        return;
    }
    
    // Copy only the upper triangle of temp to ret.
    foreach( col; 0..temp.columns ) {
        foreach( row; 0..col + 1 ) {
            ret[ row, col ] = temp[ row, col ];
        }
    } 
        
    return;
}
    
/**
Computes the Cholesky decomposition of a matrix without any copying.
Upon exiting, the upper triangle of mat will contain the result and the lower
triangle will contain zeros.  mat must be symmetric (this is not checked) but
must use general matrix storage even though it must be symmetric, since the 
lower triangle is used as scratch space.
*/
void choleskyDestructive( M )( ref M mat )
if( isMatrix!M && isGeneralMatrixStorage!( M.Storage, M.ElementType ) ) {
    assert( mat.rows == mat.columns, format(
        "mat must be symmetric for choleskyDestructive.  Got a matrix with " ~
        "%s rows, %s columns.",  mat.rows, mat.columns
    ));

    int info;    
    lapack.potrf!( 'U' )
        ( mat.length, mat.data, mat.rows, info );
        
    // TODO:  Establish a consistent standard for error handling.
    // For now, just return a matrix of NaNs.
    if( info > 0 ) {
        mat.data[ 0..mat.length ] = M.ElementType.init;
        return;
    }
    
    foreach( col; 0..mat.columns ) {
        foreach( row; col + 1..mat.columns ) {
            mat[ row, col ] = 0;
        }
    }
}
    
unittest {
    import std.stdio;
    import std.math;
    
    // Easy way to get a positive definite matrix to test with:
    auto mat = SymmetricMatrix!double( [
        [14.0, 20, 31], 
        [20.0, 62, 79], 
        [31.0, 79, 122]
    ] );
    
    auto ch = cholesky( mat );
    alias approxEqual ae;  // Save typing.
    assert( ae( ch[ 0, 0 ], 3.74166 ) );
    assert( ae( ch[ 0, 1 ], 5.34522 ) );
    assert( ae( ch[ 0, 2 ], 8.28510 ) );
    assert( ae( ch[ 1, 1 ], 5.78174 ) );
    assert( ae( ch[ 1, 2 ], 6.00412 ) );
    assert( ae( ch[ 2, 2 ], 4.16025 ) );
    
    auto alloc = newRegionAllocator();
    auto ch2 = cholesky( mat, alloc );
    foreach( row; 0..3 ) foreach( col; 0..3 ) {
        assert( ch[ row, col ] == ch2[ row, col ]);
    }
    
    Matrix!double gen;
    gen[] = mat;
    choleskyDestructive(gen);
    foreach( row; 0..3 ) foreach( col; 0..3 ) {
        assert( ch[ row, col ] == gen[ row, col ]);
    }
}
