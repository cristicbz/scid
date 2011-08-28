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

Cholesky decomposition solves for L or U.

mat may be a symmetric, triangular or general matrix, even though strictly
speaking Cholesky factorization is only defined for symmetric/Hermitian
matrices.  The triangle used to compute the Cholesky
decomposition is the upper if tri == MatrixTriangle.Upper or the lower
if tri == MatrixTriangle.Lower.  The other triangle is ignored.

Returns:  

The matrix U that would be multiplied by its conjugate transpose as
shown above to form mat, or the matrix L that would be multiplied by
its conjugate transpose to form mat, depending on the parameter tri.

Example:
---
auto mat = Matrix!double([
    [3, 1, 4],
    [5, 9, 2],
    [6, 5, 3]]
);

auto result = cholesky!(MatrixTriangle.Upper)(mat);

// The lower triangle of mat is ignored.  The above is equivalent to:
auto mat2 = SymmetricMatrix!double([
    [3, 1, 4],
    [1, 9, 2],
    [4, 2, 3]]
);

auto result2 = cholesky!(MatrixTriangle.Upper)(mat2);
---
*/
auto cholesky( MatrixTriangle tri = MatrixTriangle.Upper, M )( M mat )
if( isMatrix!M ) {  
    alias M.ElementType E;
    auto ret = TriangularMatrix!( E, tri )( mat.rows, null );
    choleskyImpl( mat, ret );
    return ret;
}

/// Ditto
auto cholesky( MatrixTriangle tri = MatrixTriangle.Upper, M, Allocator )
( M mat, Allocator alloc ) 
if( isMatrix!M && isAllocator!Allocator ) {
    alias M.ElementType E;
    alias TriangularMatrix!( E, tri ).Temporary Temp;
    auto ret = Temp( mat.rows, alloc );
    choleskyImpl( mat, ret );
    return ret;
}
     
private void choleskyImpl( M1, M2 )( ref M1 mat, ref M2 ret ) {
    alias M1.ElementType E;
    
    auto alloc = newRegionAllocator();
    auto temp = ExternalMatrixView!E( mat.rows, mat.rows, alloc );
    
    static void copyUpperTriangle( M1, M2 )( ref M1 src, ref M2 dest ) {
        foreach( col; 0..src.columns ) {
            foreach( row; 0..col + 1 ) {
                dest[ row, col ] = src[ row, col ];
            }
        }        
    }
    
    static void copyLowerTriangle( M1, M2 )( ref M1 src, ref M2 dest ) {
        foreach( col; 0..src.columns ) {
            foreach( row; col..src.columns ) {
                dest[ row, col ] = src[ row, col ];
            }
        }        
    }
    
    static if( ret.triangle == MatrixTriangle.Upper ) {
        copyUpperTriangle( mat, temp );
        enum uplo = 'U';
    } else {
        copyLowerTriangle( mat, temp );
        enum uplo = 'L';
    }
    
    int info;    
    lapack.potrf!( uplo )
        ( temp.length, temp.data, temp.rows, info );
    
    // TODO:  Establish a consistent standard for error handling.
    // For now, just return a matrix of NaNs.
    if( info > 0 ) {
        ret.data[ 0..ret.length ] = E.init;
        return;
    }
    
    static if( uplo == 'U' ) {
        copyUpperTriangle( temp, ret );
    } else {
        copyLowerTriangle( temp, ret );
    }
        
    return;
}
    
/**
Computes the Cholesky decomposition of a matrix without any copying.
Upon exiting, the upper or lower triangle of mat (depending on the value
of tri) will contain the result and the other triangle will contain zeros.  
mat is interpreted as a symmetric matrix and the upper or lower triangle
(depending on the value of tri) is used to compute the Cholesky decomposition
and the other triangle is ignored, similarly to cholesky().  However, mat must 
use general matrix storage, since one of the two triangles is used as scratch 
space.
*/
void choleskyDestructive
( MatrixTriangle tri = MatrixTriangle.Upper, M )( ref M mat )
if( isMatrix!M && isGeneralMatrixStorage!( M.Storage, M.ElementType ) ) {
    assert( mat.rows == mat.columns, format(
        "mat must be symmetric for choleskyDestructive.  Got a matrix with " ~
        "%s rows, %s columns.",  mat.rows, mat.columns
    ));

    int info;    
    lapack.potrf!( ( tri == MatrixTriangle.Upper ) ? 'U' : 'L' )
        ( mat.length, mat.data, mat.leading, info );
        
    // TODO:  Establish a consistent standard for error handling.
    // For now, just return a matrix of NaNs.
    if( info > 0 ) {
        mat.data[ 0..mat.length ] = M.ElementType.init;
        return;
    }
    
    foreach( col; 0..mat.columns ) {

        static if( tri == MatrixTriangle.Upper ) {
            immutable start = col + 1;
            immutable end = mat.columns;
        } else {
            immutable start = 0;
            immutable end = col;
        }
        
        foreach( row; start..end ) {
            mat[ row, col ] = 0;
        }
    }
}
    
unittest {
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
    
    auto ch3 = cholesky!( MatrixTriangle.Lower )( mat );
    foreach( row; 0..3 ) foreach( col; 0..3 ) {
        assert( ch3[ col, row ] == ch2[ row, col ]);
    }
    
    Matrix!double gen;
    gen[] = mat;
    choleskyDestructive( gen );
    foreach( row; 0..3 ) foreach( col; 0..3 ) {
        assert( ch[ row, col ] == gen[ row, col ]);
    }
    
    Matrix!double gen2;
    gen2[] = mat;
    choleskyDestructive!( MatrixTriangle.Lower )( gen2 );
    foreach( row; 0..3 ) foreach( col; 0..3 ) {
        assert( gen2[ col, row ] == gen[ row, col ]);
    }
}
