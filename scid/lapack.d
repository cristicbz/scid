module scid.lapack;
import scid.blas;
import scid.common.traits, scid.common.meta;
import std.algorithm, std.math, std.conv;
import std.ascii, std.exception;

// debug = lapackCalls;
//version = nodeps;

debug( lapackCalls ) {
	import std.stdio;
	import scid.internal.assertmessages;
}



version( nodeps ) {
	private enum forceNaive = true;
} else {
	private enum forceNaive = false;
	static import scid.bindings.lapack.dlapack;
	alias scid.bindings.lapack.dlapack lapack_;

}

// Save typing
int toi(size_t x) { return to!int(x); }

struct lapack {
	static void laswp( T )( size_t n, T *a, size_t lda, size_t k1, size_t k2, int *ipiv, size_t incx ) {
		debug( lapackCalls )
			writef( "laswp( %s, %s, %s, %s ) ", matrixToString( 'n', n, n, a, lda ), k1, k2, stridedToString( ipiv, n, incx ) );

		naive_.laswp( n, a, lda, k1, k2, ipiv, incx );

		debug( lapackCalls ) {
			writeln( "=> ", matrixToString( 'n', n, n, a, lda ) );
		}

	}

	static void getrs( char trans, T )( size_t n, size_t nrhs, T *a, size_t lda, int *ipiv, T *b, size_t ldb, ref int info ) {
		debug( lapackCalls )
			writef( "getrs( %s, %s, %s ) ", matrixToString( trans, n, n, a, lda ), ipiv[ 0 .. n ], matrixToString( trans, n, nrhs, b, ldb ) );

		static if( isFortranType!T && !forceNaive )
			lapack_.getrs( trans, toi(n), toi(nrhs), a, toi(lda), ipiv, b, toi(ldb), info );
		else
			naive_.xgetrs!(trans, 'L')( n, nrhs, a, lda, ipiv, b, ldb, info );

		debug( lapackCalls )
			writeln( "=> ", matrixToString( trans, n, nrhs, b, ldb ) );
	}

	static void gesv( T )( size_t n, size_t nrhs, T *a, size_t lda, int* ipiv, T *b, size_t ldb, ref int info ) {
		debug( lapackCalls )
			writef( "gesv( %s, %s, %s ) ", matrixToString( trans, n, n, a, lda ), matrixToString( trans, n, nrhs, b, ldb ) );

		static if( isFortranType!T && !forceNaive )
			lapack_.gesv( toi(n), toi(nrhs), a, toi(lda), ipiv, b, toi(ldb), info );
		else
			naive_.gesv( n, nrhs, a, lda, ipiv, b, ldb, info );

		debug( lapackCalls )
			writeln( "=> ", matrixToString( trans, n, nrhs, b, ldb ) );
	}

	static void getrf( T )( size_t m, size_t n, T* a, size_t lda, int* ipiv, ref int info ) {
		debug( lapackCalls )
			writef( "getrf( %s ) ", matrixToString( 'N', m, n, a, lda ) );

		static if( isFortranType!T && !forceNaive )
			lapack_.getrf( toi(m), toi(n), a, toi(lda), ipiv, info );
		else
			naive_.getrf( m, n, a, lda, ipiv, info );

		debug( lapackCalls )
			writeln( "=> ", matrixToString( 'N', m, n, a, lda ) );
	}

	static void trtri( char uplo, char diag, T )( size_t n, T *a, size_t lda, ref int info ) {
		debug( lapackCalls )
			writef( "trtri( %s, %s ) ", uplo, matrixToString( 'N', n, n, a, lda ) );

		static if( isFortranType!T && !forceNaive )
			lapack_.trtri( uplo, diag, toi(n), a, toi(lda), info );
		else
			naive_.trtri!( uplo, diag )( n, a, lda, info );

		debug( lapackCalls )
			writeln( "=> ", matrixToString( 'N', n, n, a, lda ) );
	}

	static void getri( T )( size_t n, T* a, size_t lda, int* ipiv, T* work, size_t lwork, ref int info ) {
		debug( lapackCalls )
			writef( "getri( %s ) ", matrixToString( 'N', n, n, a, lda ) );

		static if( isFortranType!T && !forceNaive )
			lapack_.getri( toi(n), a, toi(lda), ipiv, work, toi(lwork), info );
		else
			naive_.getri( n, a, lda, ipiv, work, lwork, info );

		debug( lapackCalls )
			writeln( "=> ", matrixToString( 'N', n, n, a, lda ) );
	}

	static void potrf( char uplo, T )( size_t n, T* a, size_t lda, ref int info ) {
	    debug( lapackCalls )
            writef( "potrf( %s, %s, %s, %s ) ", uplo_, matrixToString( uplo, n, n, a, lda ), lda, info );

        static if( isFortranType!T && !forceNaive )
			lapack_.potrf( uplo, toi(n), a, toi(lda), info );
		else
			naive_.potrf!( uplo )( toi(n), a, toi(lda), info );
	}

	// Extended LAPACK
	static void xgetrs( char trans, char side, T )( size_t n, size_t nrhs, T *a, size_t lda, int *ipiv, T *b, size_t ldb, ref int info ) {
		debug( lapackCalls )
			writef( "xgetrs( %s, %s, %s, %s ) ", side, matrixToString( trans, n, n, a, lda ), ipiv[ 0 .. n ], matrixToString( trans, n, nrhs, b, ldb ) );

		naive_.xgetrs!( trans, side )( n, nrhs, a, lda, ipiv, b, ldb, info );

		debug( lapackCalls )
			writeln( "=> ", matrixToString( trans, n, nrhs, b, ldb ) );
	}
}

private struct naive_ {
	private static void reportNaive_() {
		debug( lapackCalls )
			write( "<n> " );
	}

	private static void reportNaiveln_() {
		debug( lapackCalls )
			writeln( "<n> ..." );
	}

	static void laswp( T )( size_t n, T *a, size_t lda, size_t k1, size_t k2, int *ipiv, sizediff_t incx ) {
		reportNaiveln_();

		// convert FORTRAN indices
		--k1; --k2;

		enforce( n >= 0 );
		enforce( incx != 0 );
		enforce( a );
		enforce( ipiv );

		if( n == 0 )
		    return;

		if( incx > 0 ) {
		    for( auto i = k1; i <= k2 ; ++ i ) {
		        int pivot = ipiv[ i ] - 1; // convert FORTRAN index
		        if( pivot != i )
		            blas.swap( n, a + i, lda, a + pivot, lda );
		    }
		} else {
		    for( sizediff_t i = k2; i >= cast(sizediff_t) k1 ; -- i ) {
		        int pivot = ipiv[ i ] - 1; // convert FORTRAN index
		        if( pivot != i )
		            blas.swap( n, a + i, lda, a + pivot, lda );
		    }
		}
	}

	static void xgetrs( char trans_, char side_, T )( size_t n, size_t nrhs, T *a, size_t lda, int *ipiv, T *b, size_t ldb, ref int info ) {
		enum trans = cast(char) toUpper( trans_ );
		enum side = cast(char) toUpper( side_ );

		reportNaiveln_();

		enforce( n >= 0 );
		enforce( nrhs >= 0 );
		enforce( lda >= max( 1, n ) );
		enforce( ldb >= max( 1, n ) );

		if( n == 0 || nrhs == 0 )
			return;

		static if( trans == 'N' ) {
			static if( side == 'R' ) {
				blas.trsm!( 'R', 'U', 'N', 'N' )( n, nrhs, One!T, a, lda, b, ldb );
				blas.trsm!( 'R', 'L', 'N', 'U' )( n, nrhs, One!T, a, lda, b, ldb );
				for( sizediff_t i = n - 1; i >= 0 ; -- i ) {
					int pivot = ipiv[ i ] - 1;
					if( pivot != i )
						blas.swap( n, b + i * ldb, 1, b + pivot * ldb, 1 );
				}
			} else if( side == 'L' ) {
				lapack.laswp( nrhs, b, ldb, 1, n, ipiv, 1 );
				blas.trsm!( 'L', 'L', 'N', 'U' )( n, nrhs, One!T, a, lda, b, ldb );
				blas.trsm!( 'L', 'U', 'N', 'N' )( n, nrhs, One!T, a, lda, b, ldb );

			}
		} else {
			static if( side == 'R' ) {
				for( auto i = 0; i < n ; ++ i ) {
					int pivot = ipiv[ i ] - 1;
					if( pivot != i )
						blas.swap( n, b + i * ldb, 1, b + pivot * ldb, 1 );
				}
				blas.trsm!( 'R', 'L', trans, 'U' )( n, nrhs, One!T, a, lda, b, ldb );
				blas.trsm!( 'R', 'U', trans, 'N' )( n, nrhs, One!T, a, lda, b, ldb );
			} else static if( side == 'L' ) {
				blas.trsm!( side, 'U', trans, 'N' )( n, nrhs, One!T, a, lda, b, ldb );
				blas.trsm!( side, 'L', trans, 'U' )( n, nrhs, One!T, a, lda, b, ldb );
				lapack.laswp( nrhs, b, ldb, 1, n, ipiv, -1 );
			}
		}
	}

	static void gesv( T )( size_t n, size_t nrhs, T *a, size_t lda, int* ipiv, T *b, size_t ldb, ref int info ) {
		reportNaiveln_();

		enforce( n >= 0 );
		enforce( nrhs >= 0 );
		enforce( lda >= max( 1, n ) );
		enforce( ldb >= max( 1, n ) );

		lapack.getrf( toi(n),toi(n), a, toi(lda), ipiv, info );
		if( info == 0 )
			lapack.getrs!'N'( toi(n), toi(nrhs), a, toi(lda), ipiv, b, toi(ldb), info );
	}

	static void getri( T )( size_t n, T* a, size_t lda, int* ipiv, T* work, size_t lwork, ref int info ) {
		reportNaiveln_();
		info = 0;

		work[ 0 ] = n * 2;
		if( lwork == -1 ) {
			info = 0;
			return;
		}

		if( n == 0 )
			return;

		T get( size_t i, size_t j ) {
			return a[ j * lda + i ];
		}

		void set( string op = "" )( T x, size_t i, size_t j ) {
			mixin("a[ j * lda + i ] "~op~"= x;");
		}

		lapack.trtri!( 'U', 'N' )( n, a, lda, info );
		if( info > 0 )
			return ;

		for( int j = toi(n) - 1; j >= 0 ; -- j ) {
			foreach( i ; j + 1 .. n ) {
				work[ i ] = get( i, j );
				set( Zero!T, i, j );
			}


			if( j < n-1 ) {
				blas.gemv!'N'( n, n - j - 1, MinusOne!T, a + (j+1)*lda,
					lda, work + j + 1, 1, One!T, a + j * lda, 1 );
			}
		}

		for( int j = toi(n) - 1 ; j >= 0 ; -- j ) {
			int pivot = ipiv[ j ] - 1; // convert from FORTRAN index
			if( pivot != j )
				blas.swap( n, a + j * lda, 1, a + pivot * lda, 1 );
		}

	}

	static void trtri( char uplo_, char diag_, T )( size_t n, T *a, size_t lda, ref int info ) {
		reportNaiveln_();

		enum char uplo = toUpper(uplo_);
		enum char diag = toUpper(diag_);

		T get( size_t i, size_t j ) {
			return a[ j * lda + i ];
		}

		void set( string op = "" )( T x, size_t i, size_t j ) {
			mixin("a[ j * lda + i ] "~op~"= x;");
		}

		T ajj;
		if( uplo == 'U' ) {
			for( size_t j = 0; j < n; j++ ) {
				if( diag == 'N' ) {
					// assert( get( j, j ) != Zero!T, "fbti: Singular matrix in inverse." );
					if( get( j, j ) == Zero!T ) {
						info = j;
						return;
					}

					set( One!T / get(j,j), j, j );
					ajj = -get( j, j );
				} else {
					ajj = MinusOne!T;
				}

				blas.trmv!( 'U', 'N', diag )( j, a, lda, a + j*lda, 1 );
				blas.scal( j, ajj, a + j*lda, 1 );
			}
		} else {
			for( sizediff_t j = n - 1 ; j >= 0 ; -- j ) {
				if( diag == 'N' ) {
					// assert( get( j, j ) != Zero!T, "fbti: Singular matrix in inverse." );
					set( One!T / get(j,j), j, j );
					ajj = -get( j, j );
				} else {
					ajj = MinusOne!T;
				}

				if( j < n-1 ) {
					blas.trmv!( 'L', 'N', diag )( n - 1 - j, a + (j+1)*lda + j+1, lda, a + j*lda + j+1, 1 );
					blas.scal( n - j, ajj, a + (j+1) * lda + j, 1 );
				}
			}
		}
	}

	static void getrf( T )( size_t m, size_t n, T* a, size_t lda, int* pivot, ref int info ) {
		reportNaiveln_();

		n = min( m, n );

		T get( size_t i, size_t j ) {
			return a[ j * lda + i ];
		}

		void set( string op = "" )( T x, size_t i, size_t j ) {
			mixin("a[ j * lda + i ] "~op~"= x;");
		}

		for( size_t k = 0; k < n; k++ ) {
			pivot[ k ] = k;
			T maxSoFar = abs( get( k, k ) );
			for( size_t j = k + 1; j < n; j++ ) {
				T cur = abs( get(j, k) );
				if( maxSoFar <= cur ) {
					maxSoFar = cur;
					pivot[ k ] = j;
				}
			}

			if( pivot[ k ] != k ) {
				foreach( j ; 0 .. n ) {
					T aux = get(k, j);
					set( get(pivot[k], j), k, j );
					set( aux, pivot[k], j );
				}
			}

			if( get(k,k) != Zero!T ) {

				foreach( i ; k + 1 .. n )
					set!"/"( get(k,k), i, k );

				foreach( i ; k + 1 .. n ) {
					foreach( j ; k + 1 .. n ) {
						set!"-"( get(i,k) * get(k,j), i, j );
					}
				}
			} else if( info == 0 ) {
				info = k + 1;
			}
			++ pivot[ k ]; // convert to FORTRAN index
		}
	}

	static void potrf( char uplo_, T )( size_t n, T* a, size_t lda, ref int info ) {
	    reportNaiveln_();import std.stdio;

        // Borrowed from Don Clugston's MathExtra library, which he
        // gave me permission to relicense under Boost.

        // This algorithm is designed for row major matrices.  Depending
        // on whether we're using the lower or upper column, transpose
        // either before or after doing the decomposition to avoid an
        // impedance mismatch.

        static if( uplo_ == 'U' ) {
            transposeSquare( a, n, lda );
        }

        ref T get( size_t i, size_t j ) nothrow {
            // This is correct because we're treating the matrix as row-
            // major.
			return a[ i * lda + j ];
		}

        foreach( i; 0..n ) {
            T sum = get(i, i);

            for( sizediff_t k = i - 1; k >= 0; --k ) {
                immutable ik = get( i, k );
                sum -= ik * ik;
            }

            auto arr1 = a[ i * lda..i * lda + i ];

            if (sum > 0.0) {
                get(i, i) = sqrt( sum );

                foreach( j; i + 1..n ) {
                    import std.numeric;
                    auto arr2 = a[ j * lda..j * lda + i ];
                    immutable dot = dotProduct( arr1, arr2 );

                    sum = get( i, j ) - dot;
                    get( j, i ) = sum / get( i, i );
                }
            } else {
                info = toi( i );
                // not positive definite (could be caused by rounding errors)
                get( i, i ) = 0;
                // make this whole row zero so they have no further effect
                foreach( j; i + 1..n ) get( j, i ) = 0;
            }
        }

        static if( uplo_ == 'L' ) {
            transposeSquare( a, n, lda );
        }
	}

	unittest {
	    import scid.matrix;
	    import std.stdio;

        auto mat = Matrix!double( [
            [14.0, 20, 31],
            [0.0, 62, 79],
            [0.0, 0, 122],
            [-1.0, -1.0, -1.0]  // Dummy row to test slicing.
        ] );

        auto sliced = mat.view( 0, 3, 0, 3 );

        int info;
        potrf!( 'U' )( 3, sliced.data, 4, info );
        assert( info == 0 );

        import std.math;
        alias approxEqual ae;

        assert( ae( sliced[ 0, 0 ], 3.74166 ) );
        assert( ae( sliced[ 0, 1 ], 5.34522 ) );
        assert( ae( sliced[ 0, 2 ], 8.28510 ) );
        assert( ae( sliced[ 1, 1 ], 5.78174 ) );
        assert( ae( sliced[ 1, 2 ], 6.00412 ) );
        assert( ae( sliced[ 2, 2 ], 4.16025 ) );
        stderr.writeln(sliced.pretty);
        auto mat2 = Matrix!double( [
            [14.0, 0, 0],
            [20.0, 62, 0],
            [31.0, 79, 122],
            [-1.0, -1.0, -1.0]  // Dummy row to test slicing.
        ] );

        auto sliced2 = mat2.view( 0, 3, 0, 3 );
        potrf!( 'L' )( 3, sliced2.data, 4, info );
        assert( info == 0 );
        foreach( i; 0..3 ) foreach( j; 0..3 ) {
            assert( ae( sliced[ i, j ], sliced2[ j, i ] ) );
        }
	}
}

// Transposes a square general matrix in physical memory.  This is useful for
// making some matrix factorizations faster by improving memory locality,
// and often has negligible overhead since most matrix factorizations are
// O(N^3) and this is O(N^2) where N is the number of rows/columns.
private void transposeSquare( T )( T* ptr, size_t n, size_t lda ) pure nothrow {

    ref T get( size_t i, size_t j ) nothrow {
        return ptr[ j * lda + i ];
    }

    foreach( i; 1..n ) {
        foreach( j; 0..i ) {
            swap( get( i, j ), get( j, i ) );
        }
    }
}

unittest {
    auto mat = [1.0, 2, 3, 4, 5, 6, 7, 8, 9];
    transposeSquare( mat.ptr, 3, 3 );
    assert( mat == [1.0, 4, 7, 2, 5, 8, 3, 6, 9] );
}

