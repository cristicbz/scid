/** Common functions used to compute matrix and vector operations.

    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/

module scid.ops.common;

public import scid.ops.eval;

import std.complex;
import std.math;
import std.traits;
import scid.ops.expression;
import scid.vector;
import scid.matrix;

import std.typecons;
import std.exception;

import scid.common.traits;

import scid.blas;
import scid.lapack;

/** Generic complex conjugate that works on both builtin complex numbers and on std.Complex */
auto gconj( T )( T z ) if( is( T : cdouble ) || is( T : cfloat ) ) {
	return conj( z );
}
/// ditto
auto gconj( T )( T z ) if( is( T : Complex!double ) || is( T : Complex!float ) ) {
	return z.conj;
}

/** Inverses the value of a Transpose value. */
template transNot( Transpose tr )             { enum transNot = cast(Transpose) !tr; }

/** Performs exclusive or on Transpose values. */
template transXor( Transpose a, Transpose b ) { enum transXor = cast(Transpose) (a ^ b); }

/** Transposes a closure (only affects RowVector <-> ColumnVector) with a given Transpose value. */
template transposeClosure( Closure A, Transpose transposed = Transpose.yes ) {
	static if( transposed ) {
		static if( A == Closure.RowVector )         enum transposeClosure = Closure.ColumnVector;
		else static if( A == Closure.ColumnVector ) enum transposeClosure = Closure.RowVector;
		else                                        enum transposeClosure = A;
	} else
		enum transposeClosure = A;
}

/** Transposes a vector type (Row <-> Column) based on the given Transpose value. */
template transposeVectorType( VectorType v, Transpose transposed = Transpose.yes ) {
	static if( transposed ) {
		static if( v == VectorType.Row )         enum transposeVectorType = VectorType.Column;
		else static if( v == VectorType.Column ) enum transposeVectorType = VectorType.Row;
		else                                     enum transposeVectorType = v;
	} else
		enum transposeVectorType = v;
}

/** Transposes a storage order (RowMajor <-> ColumnMajor) based on the given Transpose value. */
template transposeStorageOrder( StorageOrder S, Transpose transposed = Transpose.yes ) {
	static if( transposed ) {
		static if( S == StorageOrder.RowMajor )         enum transposeStorageOrder = StorageOrder.ColumnMajor;
		else static if( S == StorageOrder.ColumnMajor ) enum transposeStorageOrder = StorageOrder.RowMajor;
		else                                            enum transposeStorageOrder = S;
	} else
		enum transposeStorageOrder = S;
}

/** Transposes a vector if the given parameter is Transpose.yes. This simply means conjugate its elements if the
    element type is complex.
*/
void vectorTranspose( Transpose trans, V )( ref V vec ) {
	static if( trans && isComplexScalar!(BaseElementType!V) ) {
		auto n = vec.length;
		foreach( i ; 0 .. n )
			vec[ i ] = gconj( vec[ i ] );
	}
}

/** Mixin template that provides scaledAddition(), scale() and dot() specializations for strided vector storage
    types.
*/
mixin template StridedScalingAdditionDot() {
	import scid.common.traits : isStridedVectorStorage;
	void scaledAddition( Transpose tr = Transpose.no, S )( ElementType alpha, auto ref S rhs ) if( isStridedVectorStorage!(S, ElementType) ) {
		stridedScaledAddition!tr( alpha, rhs, this );
	}

	void scale( ElementType alpha ) {
		stridedScaling( alpha, this );
	}

	auto dot( Transpose tr = Transpose.no, Other )( auto ref Other other ) if( isStridedVectorStorage!(Other, ElementType) ) {
		return stridedDot!tr( this, other );
	}
}

/** Mixin template that provides scaledAddition() and scale() for general matrix storage types. */
mixin template GeneralMatrixScalingAndAddition() {
	private import scid.lapack;
	private import scid.internal.regionallocator;

	void scaledAddition( Transpose tr = Transpose.no, S )( ElementType alpha, auto ref S rhs ) if( isGeneralMatrixStorage!(S, ElementType) ) {
		assert( rhs.rows == this.rows && rhs.columns == this.columns, "Matrix size mismatch in addition." );
		static assert( isGeneralMatrixStorage!(typeof(this), ElementType) );
		generalMatrixScaledAddition!tr( alpha, rhs, this );
	}

	void scale( ElementType alpha ) {
		static assert( isGeneralMatrixStorage!(typeof(this), ElementType) );
		generalMatrixScaling!storageOrder( this.rows, this.columns, alpha, this.data, this.leading );
	}

	void invert() {
		import scid.lapack;
		import std.exception;

		if( this.empty )
			return;

		int info;
		size_t n = this.rows;
		assert( n == this.columns, "Inversion of non-square matrix." );

		auto alloc = newRegionAllocator();
		auto ipiv  = alloc.uninitializedArray!( int[] )( n );
		auto work  = alloc.uninitializedArray!( ElementType[] )( n );

		lapack.getrf( n, n, this.data, this.leading, ipiv.ptr, info );
		lapack.getri( n, this.data, this.leading, ipiv.ptr, work.ptr , work.length, info );

		enforce( info == 0, "Inversion of singular matrix." );
	}

	void solveRight( Transpose transM = Transpose.no, Side side, Dest )( auto ref Dest dest ) if( isStridedVectorStorage!Dest || isGeneralMatrixStorage!Dest ) {
		import scid.blas;
		import scid.lapack;
		import std.exception;

		if( this.empty && dest.empty )
			return;

		size_t n = this.rows;                                           // check that the matrix is square
		assert( n == this.columns, "Inversion of non-square matrix." );

		enum vectorRhs = isStridedVectorStorage!Dest;                   // if we're solving for a vector
		auto alloc = newRegionAllocator();                              // allocator for all the copies we need to perform
		auto ipiv = alloc.uninitializedArray!(int[])( n ).ptr;          // array to store the permutation matrix of the LU decomposition
		auto a = alloc.uninitializedArray!(ElementType[])( n * n ).ptr; // array to copy this matrix's data
		int info = 0;                                                   // error number of the LAPACK calls

		// If the storage orders are different and we need to perform a hermitian transMpose then the result is
		// the same is taking the conjugate without transMposing.
		static if( !vectorRhs )
			enum bool conjugateNoTrans = isComplexScalar!ElementType && transM && storageOrder != storageOrderOf!Dest;
		else
			enum bool conjugateNoTrans = false;

		// copy this matrix's data (don't transpose, we'll use getrs's transpose for performance)
		static if( conjugateNoTrans )
			blas.xgecopyc!'n'( n, n, this.cdata, this.leading, a, n );
		else
			blas.xgecopy!'n'( n, n, this.cdata, this.leading, a, n );

		static if( vectorRhs ) {
			ElementType *b;          // if the rhs is a vector, then make it look like a matrix, by defining
			enum nrhs = 1;           // the no. of columns (nrhs) and
			auto ldb = dest.length;   // the leading dimension

			if( dest.stride != 1 ) {
				// if the destination's stride is != 1, we need to copy it to a temp array.
				b = alloc.uninitializedArray!(ElementType[])( ldb ).ptr;
				blas.copy( ldb, dest.cdata, dest.stride, b, 1 );
			} else {
				// if the stride is 1, we can simply use dest's data for b
				b = dest.data;
			}
		} else {
			ElementType *b = dest.data; // if the rhs is a matrix simply grab the parameters from it
			auto nrhs = dest.columns;
			auto ldb = dest.leading;
		}

		// Figure out what kind of transposition we've got (transM == true means it could be hermitian)
		static if( vectorRhs ) {
			enum chTrans = transM ? (isComplexScalar!ElementType ? 'C' : 'T') : 'N';
		} else static if( transposeStorageOrder!(storageOrder, transM) != storageOrderOf!Dest ) {
			static if( !isComplexScalar!ElementType || conjugateNoTrans )
				enum chTrans = 'T';
			else
				enum chTrans = 'C';
		} else {
			enum chTrans = 'N';
		}

		enum chSide = (side == Side.Left) ? 'L' : 'R';

		lapack.getrf( n, n, a, n, ipiv, info );                    // perform LU decomposition
		enforce( info == 0, "Inversion of singular matrix." );
		lapack.xgetrs!(chTrans, chSide)( n, nrhs, a, n, ipiv, b, ldb, info ); // perform inv-mult

		static if( vectorRhs ) {
			// copy the data back to dest if needed
			if( dest.stride != 1 )
				blas.copy( ldb, b, 1, dest.data, dest.stride );
		}
	}

	void matrixProduct( Transpose transA = Transpose.no, Transpose transB = Transpose.no, A, B )
			( ElementType alpha, auto ref A a, auto ref B b, ElementType beta ) if( isGeneralMatrixStorage!A && isGeneralMatrixStorage!B ) {
		import scid.blas;

		// Are the matrices row major?
		enum aRowMajor = storageOrderOf!A == StorageOrder.RowMajor;
		enum bRowMajor = storageOrderOf!B == StorageOrder.RowMajor;
		enum cRowMajor = storageOrder == StorageOrder.RowMajor;

		// Should a and be read minor-first or major-first
		enum aByMinor = transA ^ aRowMajor ^ cRowMajor;
		enum bByMinor = transB ^ bRowMajor ^ cRowMajor;

		// Do the matrices have complex elements?
		enum complexElements = isComplexScalar!ElementType;

		static if( complexElements ) {
			enum aConj = transA && !aByMinor, aHerm = transA && aByMinor;
			enum bConj = transB && !bByMinor, bHerm = transB && bByMinor;

			static if( aByMinor )
				enum aTransChar = bConj || !transA ? 'T' : 'C';
			else
				enum aTransChar = 'N';

			static if( bByMinor )
				enum bTransChar = aConj || !transB ? 'T' : 'C';
			else
				enum bTransChar = 'N';

			enum conjugateResult = aConj && bHerm || bConj && aHerm || aConj && bConj;

			enum aTemp = aConj && !transB;
			enum bTemp = bConj && !transA;

			static if( aTemp || bTemp )
				auto alloc = newRegionAllocator();

			static if( aTemp ) {
				auto aData    = cast(ElementType*)alloc.allocate( a.rows * a.columns );
				auto aLeading = a.minor;

				blas.xgecopyc!'N'( a.minor, a.major, a.cdata, a.leading, aData, aLeading );
			} else {
				auto aData = a.cdata, aLeading = a.leading;
			}

			static if( bTemp ) {
				auto bData    = cast(ElementType*)alloc.allocate( b.rows * b.columns );
				auto bLeading = b.minor;

				blas.xgecopyc!'N'( b.minor, b.major, b.cdata, b.leading, bData, bLeading );
			} else {
				auto bData = b.cdata, bLeading = b.leading;
			}
		} else {
			enum conjugateResult = false;
			auto aData = a.cdata, aLeading = a.leading;
			auto bData = b.cdata, bLeading = b.leading;

			enum aTransChar = aByMinor ? 'T' : 'N';
			enum bTransChar = bByMinor ? 'T' : 'N';
		}

		static if( !cRowMajor ) {
			static if( !aByMinor )
				size_t m = a.minor, k = a.major;
			else
				size_t m = a.major, k = a.minor;

			static if( !bByMinor )
				size_t n = b.major;
			else
				size_t n = b.minor;

			if( !beta )
				this.resize( m, n, null );

			blas.gemm!( aTransChar, bTransChar )( m, n, k, alpha, aData, aLeading, bData, bLeading, beta, this.data, this.leading );
		} else {
			static if( !bByMinor )
				size_t m = b.minor, k = b.major;
			else
				size_t m = b.major, k = b.minor;

			static if( !aByMinor )
				size_t n = a.major;
			else
				size_t n = a.minor;

			if( !beta )
				this.resize( m, n, null );

			blas.gemm!( bTransChar, aTransChar )( m, n, k, alpha, bData, bLeading, aData, aLeading, beta, this.data, this.leading );
		}

		static if( conjugateResult )
			blas.xgecopyc!'N'(m,n,this.data,this.leading,this.data,this.leading);
	}
}

/** Compute the dot product of a row and a column of possibly transposed matrices. */
auto rowColumnDot( Transpose transA = Transpose.no, Transpose transB = Transpose.no, A, B )( auto ref A a, size_t i, auto ref B b, size_t j ) {
	static if( !transA && !transB ) {
		return evalDot!( transA, transB )( a.row( i ), b.column( j ) );
	} else static if( !transA && transB ) {
		return evalDot!( transA, transB )( a.row( i ), b.row( j ) );
	} else static if( transA && !transB ) {
		return evalDot!( transA, transB )( a.column( i ), b.column( j ) );
	} else {
		return evalDot!( transA, transB )( b.column( i ), b.row( j ) );
	}
}

/** Specialized copy for strided vector storage types. */
void stridedCopy( Transpose tr, Source, Dest )( auto ref Source source, auto ref Dest dest ) if( isStridedVectorStorage!(Source,BaseElementType!Dest) && isStridedVectorStorage!Dest ) {
	dest.resize( source.length, null );
	static if( isComplexScalar!(BaseElementType!Source) && tr ) {
		blas.xcopyc( dest.length, source.cdata, source.stride, dest.data, dest.stride );
	} else
		blas.copy( dest.length, source.cdata, source.stride, dest.data, dest.stride );

}

/** Specialized scaled addition (axpy) for strided vector storage types. */
void stridedScaledAddition( Transpose tr, Scalar, Source, Dest )( Scalar alpha, ref Source source, auto ref Dest dest ) if( isStridedVectorStorage!(Source,Scalar) && is( Scalar : BaseElementType!Dest ) ) {
	assert( source.length == dest.length, "Length mismatch in strided addition." );
	static if( isComplexScalar!(BaseElementType!Source) && tr )
		blas.xaxpyc( dest.length, alpha, source.cdata, source.stride, dest.data, dest.stride );
	else
		blas.axpy( dest.length, alpha, source.cdata, source.stride, dest.data, dest.stride );
}

/** Specialized scaling for strided vector storage types. */
void stridedScaling( Scalar, Dest )( Scalar alpha, auto ref Dest dest ) if( isStridedVectorStorage!Dest && is( Scalar : BaseElementType!Dest ) ) {
	blas.scal( dest.length, alpha, dest.data, dest.stride );
}

/** Specialized dot for strided vector storage types. */
auto stridedDot( Transpose transA, A, B )( auto ref A a, auto ref B b ) if( isStridedVectorStorage!(A,BaseElementType!B) && isStridedVectorStorage!B ) {
	assert( a.length == b.length, "Length mismatch in strided dot." );
	static if( !isComplexScalar!(BaseElementType!A) )
		auto r = blas.dot( a.length, a.cdata, a.stride, b.cdata, b.stride );
	else static if( transA )
		auto r = blas.dotc( a.length, a.cdata, a.stride, b.cdata, b.stride );
	else
		auto r = blas.dotu( a.length, a.cdata, a.stride, b.cdata, b.stride );
	return r;
}

/** Specialized scaling for general matrix storage types. */
void generalMatrixScaling( StorageOrder order, E )( size_t rows, size_t cols, E alpha,  E* data, size_t leading ) {
	static if( order == StorageOrder.RowMajor ) {
		alias rows major;
		alias cols minor;
	} else {
		alias rows minor;
		alias cols major;
	}

	if( leading == minor ) {
		// scale everything
		blas.scal( rows * cols, alpha, data, 1 );
	} else {
		// scale by-major
		auto e = data + major * leading;
		for( ; data < e ; data += leading )
			blas.scal( minor, alpha, data, 1 );
	}
}

/** Specialized copying for general matrix storage types. */
void generalMatrixCopy( Transpose tr, S, D )( auto ref S source, auto ref D dest ) {
	alias BaseElementType!S T;
	alias transposeStorageOrder!(storageOrderOf!S, tr) srcOrder;
	alias storageOrderOf!D dstOrder;

	static if( !isComplexScalar!T ) {
		blas.xgecopy!( (srcOrder == dstOrder) ? 'N' : 'T')(
			dest.rows,
			dest.columns,
			source.cdata,
			source.leading,
			dest.data,
			dest.leading
		);
	} else static if( isComplexScalar!T ) {
		static if( srcOrder == dstOrder ) {
			static if( tr ) {
				blas.xgecopyc!'N'( dest.rows, dest.columns,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			} else {
				blas.xgecopy!'N'( dest.rows, dest.columns,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			}
		} else {
			static if( tr ) {
				blas.xgecopy!'C'( dest.rows, dest.columns,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			} else {
				blas.xgecopy!'T'( dest.rows, dest.columns,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			}
		}
	}
}

/** Specialized scaled addition (axpy) for general matrix storage types. */
void generalMatrixScaledAddition( Transpose tr, S, D, E )( T alpha, auto ref S source, auto ref D dest ) {
	alias transposeStorageOrder!(storageOrderOf!S, tr) srcOrder;
	alias storageOrderOf!D dstOrder;

	static if( !isComplexScalar!T ) {
		blas.xgeaxpy( (srcOrder == dstOrder) ? 't' : 'n',
			dest.rows, dest.columns, alpha,
			source.cdata, source.leading,
			dest.data, dest.leading
		);
	} else static if( isComplexScalar!T ) {
		static if( srcOrder == dstOrder ) {
			static if( tr ) {
				blas.xgeaxpyc( 'n', dest.rows, dest.columns, alpha,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			} else {
				blas.xgeaxpy( 'n', dest.rows, dest.columns, alpha,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			}
		} else {
			static if( tr ) {
				blas.xgeaxpy( 'c', dest.rows, dest.columns, alpha,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			} else {
				blas.xgeaxpy( 't', dest.rows, dest.columns, alpha,
					source.cdata, source.leading,
					dest.data,    dest.leading
				);
			}
		}
	}
}

/**
This function takes a vector and returns an ExternalMatrixView of
the vector.  The purpose of this is to support things like vector-matrix
multiplication by treating the vector on the LHS as a matrix.  It hacks
around the ref counting system because it's only supposed to be used
internally for evaluating e.g. vector-matrix multiplication.

If V is already a matrix and not a vector, this is a no-op.
*/
auto toMatrix( V )( V vec ) {
    alias ExternalMatrixView!( typeof( vec[0] ) ) Ret;
    immutable n = vec.length;
    auto ptr = cast(Unqual!(typeof(*vec.cdata))*) vec.cdata;

    static if( closureOf!V == Closure.ColumnVector ) {
        return Ret( 1, ptr[0..n] );
    } else {
        return Ret( n, ptr[0..n] );
    }
}

unittest {
    auto col = Vector!double([1, 2, 3]);
    auto row = eval(col.t);

    auto colMat = toMatrix(col);
    auto rowMat = toMatrix(row);

    assert(colMat[0, 0] == 1);
    assert(colMat[1, 0] == 2);
    assert(colMat[2, 0] == 3);

    assert(rowMat[0, 0] == 1);
    assert(rowMat[0, 1] == 2);
    assert(rowMat[0, 2] == 3);

    assert(colMat.cdata is rowMat.cdata);
    assert(rowMat.cdata is col.cdata);
}
