module scid.demo;

version( demo ) {
	import scid.matvec;
	import scid.storage.diagonalmat;
	import scid.common.traits, scid.common.meta;
	import scid.internal.regionallocator;
	
	import std.stdio, std.conv, std.typetuple;
	import std.complex, std.range, std.exception;
	import std.string, std.math;
	
	void main() {
        auto x = DiagonalMatrix!double([1.,2,3,4]), y=DiagonalMatrix!double([2.,2.,2.,2.]);
        //  auto w = SymmetricMatrix!double([[1.,2.,3.,4],[1.,2.,3.,4],[1.,2.,3.,4],[1.,2.,3.,4]]);
        eval( x[0..2][] * y[][0..2] );
        //writeln(z.pretty);

		testIssue51();
		readln();
	}
	
	// Methods demonstrating basic functionality.
	
	/** Syntax for basic expressions. */
	void basicExpressions()() {
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
	
	/** Using vectors and matrices and ranges. */
	void rangeInterface()() {
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
	
	
	/** Directly accessing the data of vectors and matrices. */
	void dataInterface()() {
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
	
	/** Creating views of existing data that allow it to be treated as a matrix or vector. */
	void externalViews()() {
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
	
	// ====================================================== TESTING ===========================================
	// These should be moved into unittests at some point, for now OPTLINK makes it more comfortable to have them
	// this way...
	
	/** myApproxEqual that works on complex elements as well */
	bool myApproxEqual( T, U )( T a, U b ) {
		static if( is( T : cdouble ) && is( U : cdouble ) )
			return myApproxEqual( a.re, b.re ) && myApproxEqual( a.im, b.im );
		else
			return approxEqual( a, b );
	}
	
	/** Enforce the data & dimensions of a matrix. */
	void enforceMatData( M, E )( auto ref M m, size_t r, size_t c, E[] expected ) {
		//debug {
			enum epsilon = 1e-3;
			enforce( m.rows == r, format("Wrong no. of rows %d vs %d", m.rows, r) );
			enforce( m.columns == c, format("Wrong no. of rows %d vs %d", m.columns, c ) );
			auto a = m.cdata[ 0 .. expected.length ];
			auto b = to!(BaseElementType!M[])(expected);
			foreach( i ; 0 .. a.length ) {
				enforce( myApproxEqual( a[i], b[i] ),
				   "Expected " ~ to!string(expected) ~ ", got "~ to!string(a) ~
					" (" ~ to!string(i) ~  ")"
				);
			}
		//}
	}
	
	/** Enforce the data and length of a vector. */
	void enforceVecData( V, E )( auto ref V v, E[] expected ) {
		//debug {
			enforce( v.length == expected.length, format("Wrong vector length: %d vs. %d", v.length, expected.length) );
			enforce( v.cdata[ 0 .. expected.length ] == to!(BaseElementType!V)(expected) );
		//}
	}
	
	/** Enforce that a matrix is myApproxEqual to a built-in matrix. */
	void enforceMatApproxEqual( int line = __LINE__, A, E )( auto ref A a, E[][] b ) {
		string message() { return "At " ~ to!string(line) ~ ": Expected " ~ to!string(b) ~ ", got " ~ a.toString() ~ ". "; }
			
		enforce( a.rows == b.length, message() ~ "First dimensions do not match." );
			
		if( a.rows == 0 )
			return;
			
		enforce( a.columns == b[0].length, message() ~ "Second dimensions do not match." );
			
		foreach( i ; 0 .. a.rows ) {
			foreach( j ; 0 .. a.columns ) {
				enforce( myApproxEqual( a[ i, j ], b[ i ][ j ] ), message() ~
					format("Elements at position (%d,%d) do not match.", i, j)
				);
			}
		}
	}
	
	// Some of these types have to be disabled otherwise the compiler runs out of memory.
	/** A TypeTuple of all matrix types to be tested by generic tests. */
	template MatrixTypes( T ) {
		alias TypeTuple!(
			//Matrix!T,
			//Matrix!(T,StorageOrder.RowMajor),
			//TriangularMatrix!T,
			//TriangularMatrix!(T, MatrixTriangle.Lower ),
			//TriangularMatrix!(T, MatrixTriangle.Upper, StorageOrder.RowMajor),
			//TriangularMatrix!(T, MatrixTriangle.Lower, StorageOrder.RowMajor),
			//SymmetricMatrix!T,
			// SymmetricMatrix!(T, MatrixTriangle.Lower ),
			// SymmetricMatrix!(T, MatrixTriangle.Upper, StorageOrder.RowMajor),
			// SymmetricMatrix!(T, MatrixTriangle.Lower, StorageOrder.RowMajor),
			DiagonalMatrix!T
		) MatrixTypes;
	}
	
	/** Check that empty matrices and vectors don't cause crashes. */
	void emptyMatVecTest()() {
		Matrix!double m1,m2;
		Vector!double v1,v2;
		
		// We're really checking there's no assertions or memory problems
		// Empty vectors
		v1[] = v2;
		v2[] = v1[ 0 .. 0 ];
		v1[] *= 2.0;
		v1[] += v2;
		double x = eval( v1.t*v2 );
			
		// Empty matrices
		m1[ 0 .. 0 ][ 0 .. 0 ] = m2;
		m1[] += m2;
		m1[] *= m2;
		
		// The inverse of an empty matrix is the empty matrix. After all, it satisfies
		// the properties of the identity matrix.
		enforce( eval( inv( m1 ) ).empty );
		enforce( ( eval( inv( m1 ) * m2 ) ).empty );
	}
	
	/** Syntactically test all the operations on all the matrix types. */
	void syntaxTest()() {
		alias TypeTuple!(double) ElementTypes;
		foreach(T; ElementTypes) {
			enum z = One!T;
			T[][] minit = [[z, z, z], [z, z, z], [z, z, z]];
			foreach(LhsType; MatrixTypes!T) {
				auto lhs = LhsType( minit );
				foreach(RhsType; MatrixTypes!T) {
					auto rhs = RhsType( minit );
					
					eval(lhs + rhs*z);
					eval(lhs*z - rhs);
					eval(lhs * rhs*z);
					eval(lhs.column(0)*z + rhs.column(1));
					eval(lhs.row(0) - rhs.row(1)*z);
					eval(lhs.row(0) * rhs.column(0)*z);
					
					eval( lhs.t*(lhs + rhs*z) );
					eval( (lhs.column(0) - rhs.column(0)*z).t*lhs );
					eval( (lhs.column(0)*z + rhs.column(0)).t*lhs.column(1) );
					
					lhs[] += rhs;
					lhs[] -= rhs;
					lhs[] *= rhs;
					lhs[] *= z;
					lhs[] /= z;
					
					lhs[0][] *= z;
					lhs[0][] /= z;
					lhs[0][] += lhs[][0].t;
					lhs[][0] -= lhs[1][].t;
					
					lhs[] += z;
					lhs[] -= z;
					lhs[] = z;
					lhs[] = (lhs + z)*(lhs - lhs[0][]*lhs[][0]);
					
					lhs[1..3][1..3] *= rhs[0..2][0..2];
				}
			}
		}
	}
	
	/** Test a medley of matrix operations with double elements. */
	void dMatOpsTest()() {
		alias Matrix!double            dGeMat;
		alias SymmetricMatrix!double   dSyMat;
		
		auto a = dGeMat( 3, [1.,2,3,4,5,6,7,8,9] );
		auto b = dGeMat( 3, [1.,2,3,4,5,6] );
		
		dGeMat c = b * a;
		enforceMatData( c, 2, 3, [22,28,49,64,76,100] );
		c[] = c[0 .. 2][ 0 .. 2 ].t * ( (b[][0] - a[1..3][0]).t * eval(c[][0]) ) / 50.;
		enforceMatData( c, 2, 2, [-22,-49,-28,-64] );
		
		dSyMat s = c.t*c;
		
		enforceMatData( s, 2, 2, [2885, 3752, 4880] );
		
		auto d = eval( s - dSyMat([2800.,3700,4800]) );
		static assert( is( typeof(d) : dSyMat ) );
		enforceMatData( d, 2, 2, [85,52,80] );
		enforce( d[1][0] == 52 );
		
		auto e = eval( d - b[0..2][1..3]*10 );
		static assert( is( typeof(e) : dGeMat ) );
		enforceMatData( e, 2, 2, [ 55, 12, 2, 20 ] );
	}
	
	/** Test inversions on matrices of doubles. */
	void dMatInvTest()() {
		auto x = Matrix!double([ [ 1, 2, 3 ], [  1, 1, 1 ], [ 2, -1, 1 ] ]);
		auto y = Matrix!double([ [ 1, 2, 0 ], [ -2, 3, 4 ], [ 0,  2, 1 ] ]);
		
		// Thank you Octave...
		enforceMatData( eval(inv(x)*y), 3, 3,     [  -2.4,  -2.2,   2.6,   2.6,   1.8,  -1.4,   4.2,   3.6,  -3.8 ] );
		enforceMatData( eval(y*inv(x)), 3, 3,     [  -0.8,   2.6,   0.2,   3.0,  -3.0,   1.0,  -0.6,  -0.8,  -0.6 ] );
		enforceMatData( eval(inv(x.t)*y), 3, 3,   [   0.0,  -1.0,   1.0,  -0.2,   3.0,  -0.4,  -0.2,   3.0,  -1.4 ] );
		enforceMatData( eval(y*inv(x.t)), 3, 3,   [   1.6,   4.6,   2.2,   1.8,   1.8,   1.6,  -1.4,  -3.4,  -1.8 ] );
		enforceMatData( eval(inv(y)*x), 3, 3,     [  -9.0,   5.0,  -8.0,  20.0,  -9.0,  17.0,   9.0,  -3.0,   7.0 ] );
		enforceMatData( eval(x*inv(y)), 3, 3,     [  13.0,   7.0,  16.0,   6.0,   3.0,   7.0, -21.0, -11.0, -27.0 ] );
		enforceMatData( eval(inv(y.t)*x), 3, 3,   [  11.0,   5.0, -18.0,   4.0,   1.0,  -5.0,  17.0,   7.0, -27.0 ] );
		enforceMatData( eval(x*inv(y.t)), 3, 3,   [ -15.0,  -1.0,   0.0,   8.0,   1.0,   1.0, -13.0,  -1.0,  -1.0 ] );
		
		enforceMatData( eval(x.t*inv(y.t)), 3, 3, [  -9.0,  20.0,   9.0,   5.0,  -9.0,  -3.0,  -8.0,  17.0,   7.0 ] );
		enforceMatData( eval(y.t*inv(x.t)), 3, 3, [  -2.4,   2.6,   4.2,  -2.2,   1.8,   3.6,   2.6,  -1.4,  -3.8 ] );
	}
	
	/** Specifically test all cases of matrix products for doubles (32 cases...). */
	void dMatProdTest()() {
		auto aInit    = [ [ 1.,  2, 1 ], [ 3., -1, 2]  ];
		auto bInit    = [ [ 1., -2 ],    [ -3., 2 ], [ -2., 1. ] ];
		auto correct1 = [ [-7., 3], [2., -6] ];
		auto correct2 = [ [-5.,3,1],[4.,-8,-5],[-3.,1,0] ];	
		
		alias Matrix!double                            ColMat;
		alias Matrix!( double, StorageOrder.RowMajor ) RowMat;
		
		// check all cases except mixed transpose (dimensions wouldn't match)
		auto colA = ColMat( aInit ), colB = ColMat( bInit );
		auto rowA = RowMat( aInit ), rowB = RowMat( bInit );
		ColMat colC;
		RowMat rowC;
		
		colC[] = colA * colB;     enforceMatApproxEqual( colC, correct1 );
		colC[] = colA * rowB;     enforceMatApproxEqual( colC, correct1 );
		colC[] = rowA * colB;     enforceMatApproxEqual( colC, correct1 );
		colC[] = rowA * rowB;     enforceMatApproxEqual( colC, correct1 );
		
		colC[] = colA.t * colB.t; enforceMatApproxEqual( colC, correct2 );
		colC[] = colA.t * rowB.t; enforceMatApproxEqual( colC, correct2 );
		colC[] = rowA.t * colB.t; enforceMatApproxEqual( colC, correct2 );
		colC[] = rowA.t * rowB.t; enforceMatApproxEqual( colC, correct2 );
		
		rowC[] = colA * colB;     enforceMatApproxEqual( rowC, correct1 );
		rowC[] = colA * rowB;     enforceMatApproxEqual( rowC, correct1 );
		rowC[] = rowA * colB;     enforceMatApproxEqual( rowC, correct1 );
		rowC[] = rowA * rowB;     enforceMatApproxEqual( rowC, correct1 );
		
		rowC[] = colA.t * colB.t; enforceMatApproxEqual( rowC, correct2 );
		rowC[] = colA.t * rowB.t; enforceMatApproxEqual( rowC, correct2 );
		rowC[] = rowA.t * colB.t; enforceMatApproxEqual( rowC, correct2 );
		rowC[] = rowA.t * rowB.t; enforceMatApproxEqual( rowC, correct2 );
		
		// check mixed transpose
		auto correct3 = [[-8.,4],[5.,-6]];
		auto correct4 = [[-3.,1],[5.,-11]];
		colA[] = colA[0..2][0..2]; colB[] = colB[0..2][0..2];
		rowA[] = rowA[0..2][0..2]; rowB[] = rowB[0..2][0..2];
		
		colC[] = colA.t * colB; enforceMatApproxEqual( colC, correct3 );
		colC[] = colA.t * rowB; enforceMatApproxEqual( colC, correct3 );
		colC[] = rowA.t * colB; enforceMatApproxEqual( colC, correct3 );
		colC[] = rowA.t * rowB; enforceMatApproxEqual( colC, correct3 );
		
		colC[] = colA * colB.t; enforceMatApproxEqual( colC, correct4 );
		colC[] = colA * rowB.t; enforceMatApproxEqual( colC, correct4 );
		colC[] = rowA * colB.t; enforceMatApproxEqual( colC, correct4 );
		colC[] = rowA * rowB.t; enforceMatApproxEqual( colC, correct4 );
		
		rowC[] = colA.t * colB; enforceMatApproxEqual( rowC, correct3 );
		rowC[] = colA.t * rowB; enforceMatApproxEqual( rowC, correct3 );
		rowC[] = rowA.t * colB; enforceMatApproxEqual( rowC, correct3 );
		rowC[] = rowA.t * rowB; enforceMatApproxEqual( rowC, correct3 );
		
		rowC[] = colA * colB.t; enforceMatApproxEqual( rowC, correct4 );
		rowC[] = colA * rowB.t; enforceMatApproxEqual( rowC, correct4 );
		rowC[] = rowA * colB.t; enforceMatApproxEqual( rowC, correct4 );
		rowC[] = rowA * rowB.t; enforceMatApproxEqual( rowC, correct4 );
	}
	
	void zMatProdTest()() {
		alias Matrix!cdouble                            ColMat;
		alias Matrix!( cdouble, StorageOrder.RowMajor ) RowMat;
		
		auto aInit    = [ [  1. + 1.i,  2. - 1.i ], [ -1. + 2.i,  1. + 1.i ] ];
		auto bInit    = [ [ -1. - 1.i,  2. - 2.i ], [  2. + 1.i,  0. - 1.i ] ];
		auto correct1 = [ [  5. - 2.i,  3. - 2.i ], [  4. + 2.i,  3. + 5.i ] ];
		auto correct2 = [ [  2. - 4.i,  3. - 4.i ], [  1. + 1.i,  6. + 1.i ] ];
		auto correct3 = [ [ -2. - 5.i, -2. - 3.i ], [  2. - 4.i,  5. - 3.i ] ];
		auto correct4 = [ [  4. + 2.i,  4. + 3.i ], [ -1. + 1.i, -1. + 6.i ] ];
		
		auto colA = ColMat( aInit ), colB = ColMat( bInit );
		auto rowA = RowMat( aInit ), rowB = RowMat( bInit );
		ColMat colC;
		RowMat rowC;
		
		// check all cases
		colC[] = colA * colB;     enforceMatApproxEqual( colC, correct1 );
		colC[] = colA * rowB;     enforceMatApproxEqual( colC, correct1 );
		colC[] = rowA * colB;     enforceMatApproxEqual( colC, correct1 );
		colC[] = rowA * rowB;     enforceMatApproxEqual( colC, correct1 );
		
		colC[] = colA.t * colB.t; enforceMatApproxEqual( colC, correct2 );
		colC[] = colA.t * rowB.t; enforceMatApproxEqual( colC, correct2 );
		colC[] = rowA.t * colB.t; enforceMatApproxEqual( colC, correct2 );
		colC[] = rowA.t * rowB.t; enforceMatApproxEqual( colC, correct2 );
		
		rowC[] = colA * colB;     enforceMatApproxEqual( rowC, correct1 );
		rowC[] = colA * rowB;     enforceMatApproxEqual( rowC, correct1 );
		rowC[] = rowA * colB;     enforceMatApproxEqual( rowC, correct1 );
		rowC[] = rowA * rowB;     enforceMatApproxEqual( rowC, correct1 );
		
		rowC[] = colA.t * colB.t; enforceMatApproxEqual( rowC, correct2 );
		rowC[] = colA.t * rowB.t; enforceMatApproxEqual( rowC, correct2 );
		rowC[] = rowA.t * colB.t; enforceMatApproxEqual( rowC, correct2 );
		rowC[] = rowA.t * rowB.t; enforceMatApproxEqual( rowC, correct2 );
		
		colC[] = colA.t * colB;   enforceMatApproxEqual( colC, correct3 );
		colC[] = colA.t * rowB;   enforceMatApproxEqual( colC, correct3 );
		colC[] = rowA.t * colB;   enforceMatApproxEqual( colC, correct3 );
		colC[] = rowA.t * rowB;   enforceMatApproxEqual( colC, correct3 );
								  
		colC[] = colA * colB.t;   enforceMatApproxEqual( colC, correct4 );
		colC[] = colA * rowB.t;   enforceMatApproxEqual( colC, correct4 );
		colC[] = rowA * colB.t;   enforceMatApproxEqual( colC, correct4 );
		colC[] = rowA * rowB.t;   enforceMatApproxEqual( colC, correct4 );
								  
		rowC[] = colA.t * colB;   enforceMatApproxEqual( rowC, correct3 );
		rowC[] = colA.t * rowB;   enforceMatApproxEqual( rowC, correct3 );
		rowC[] = rowA.t * colB;   enforceMatApproxEqual( rowC, correct3 );
		rowC[] = rowA.t * rowB;   enforceMatApproxEqual( rowC, correct3 );
								  
		rowC[] = colA * colB.t;   enforceMatApproxEqual( rowC, correct4 );
		rowC[] = colA * rowB.t;   enforceMatApproxEqual( rowC, correct4 );
		rowC[] = rowA * colB.t;   enforceMatApproxEqual( rowC, correct4 );
		rowC[] = rowA * rowB.t;   enforceMatApproxEqual( rowC, correct4 );
	}
	
	/** Test a medley of operations on complex-valued matrices. */
	void zMatOpsTest()() {
		alias Matrix!cdouble            zGeMat;
		alias SymmetricMatrix!cdouble   zSyMat;
		
		auto a = zGeMat( 3, [1.+4.i,2+3.i,3+2.i,4+1.i,5+0.i,6-1.i,7-2.i,8-3.i,9-4.i] );
		auto b = zGeMat( 3, [1.+2.i,2.+1.i,3+0.i,4-1.i,5-2.i,6-3.i] );
				   
		zGeMat c = b * a;
		enforceMatData( c, 2, 3, [ 18 + 19.i, 33 + 22i, 45 - 8i, 60 - 23i, 72 -35i, 87 - 68i ] );
		
		
		c[] = c[0 .. 2][ 0 .. 2 ].t * ( (b[][0] - a[1..3][0]).t * eval(c[][0]) ) / (10.+0.i);
		enforceMatData( c, 2, 2, [-146.60 + 192.80i, -422.00 -  28.60i, -281.60 + 235.40i, -575.00 - 151.60i] );
		
		c[] = c + zGeMat([[150-190i,280-230i], [430+28i,570+150i]]);
		zSyMat s = c.t*c;
		
		enforceMatData( s, 2, 2, [ 83.760 +  0.000i,  
		                   -29.360 +  7.040i,
		                    59.280 +  0.000i ]);
		
		
		enforce( abs(s[1][0] - (-29.360 - 7.040i)) <= 1e-3 );
		
		auto d = eval( s - zSyMat([80.+0.i,-28,59]) );
		static assert( is( typeof(d) : zSyMat ) );
		
		enforceMatData( d, 2, 2, [ 3.76 + 0.0i, -1.36 + 7.04i, 0.28 + 0.0i ] );
		
		
		enforce( abs(d[1][0] - (-1.36 - 7.04i)) <= 1e-3 );
		
		
		auto e = eval( d + b[0..2][1..3]*(10.+0.i) );
		static assert( is( typeof(e) : zGeMat ) );
		
		enforceMatData( e, 2, 2, [ 33.760 +  0.000i, 38.640 - 17.040i, 48.640 - 12.960i, 60.280 - 30.000i] );
	}
	
	/** Issue 52 - Inverting Scalars */
	void testIssue52()() {
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
	
	/** Issue 51 - Invalid multiplication between ColumnVector and Matrix */
	void testIssue51()() {
	    // Testing row vectors, too.
	    auto col = vector( [1.0, 2, 3] );
	    auto row = vector( [4.0, 5, 6] );
	    auto mat = matrix( [[4.0, 5, 6]] );
	    
	    auto rowRes = eval( col * row.t );
	    auto matRes = eval( col * mat );
	    
	    foreach( res; [rowRes, matRes] ) {
	        assert(	res[0, 0] == 4 );
	        assert(	res[0, 1] == 5 );
	        assert(	res[0, 2] == 6 );
	        assert(	res[1, 0] == 8 );
	        assert(	res[1, 1] == 10 );
	        assert(	res[1, 2] == 12 );
	        assert(	res[2, 0] == 12 );
	        assert(	res[2, 1] == 15 );
	        assert(	res[2, 2] == 18 );
	    }
	}
	
	/** Issue 50 - Matrix slice-slice-assign ends up transposed	*/
	void testIssue50()() {
		auto a = Matrix!double([[1.0, 2], [4.0, 5]]);
		auto b = Matrix!double(3, 3);
		b[0..2][0..2] = a;
		enforce( b[0,0] == 1. && b[0,1] == 2. && b[1,0] == 4. && b[1,1] == 5. );
	}
	
	/** Issue 49 - Wrong index-slice-assign */
	void testIssue49()() {
		auto mat = Matrix!double([[1,2,3],[4,5,6],[7,8,9]]);
		auto vec = Vector!double([100, 200]).t;
		mat[2][0..2] = vec;
		enforceMatData( mat, 3, 3, [1.,4,100,2,5,200,3,6,9] );
	}
	
	/** Issue 35 - Assignment should return value for chaining */
	void testIssue35()() {
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
	
	/** Issue 32 - D'tor problems/assert failures, Linux Only */
	void testIssue32()() {
		auto xTx = Matrix!double(1, 1);
		auto xTy = Matrix!double(1, 2);
		xTx[0, 0] = 31;

		xTy[0, 0] = 41;	
		xTy[0, 1] = 59;

		auto ret = eval(inv(xTx) * xTy);
	}
	
	/** Issue 48 - Vector should be a random-access range */
	void testIssue48()() {
		import std.range;
		static assert(isInputRange!(Vector!double));         
		static assert(isForwardRange!(Vector!double));       
		static assert(isBidirectionalRange!(Vector!double)); 
		static assert(isRandomAccessRange!(Vector!double)); 
	}
	
	/** Isuue 54 - Matrix slices have copy semantics, while row/column slices have view semantics */
	void testIssue54()() {
		auto mat   = Matrix!double([ [1.,2.], [3., 4.] ]);	
		auto row   = mat[ 0 ][ 0 .. 2 ];
		auto col   = mat[ 0 .. 2 ][ 0 ];
		auto slice = mat[ 0 .. 1 ][ 0 .. 1 ];
		
		// Copy semantics of matrix slices:
		slice[ 0, 0 ] = 10.0;
		enforce( mat[ 0, 0 ] = 1. );
		
		// View semantics of row and column slices:
		row[ 0 ] = 100.;
		enforce( mat[ 0, 0 ] == 1. && col[ 0 ] == 1. );
		
		col[ 0 ] = 200.;
		enforce( mat[ 0, 0 ] == 1. && row[ 0 ] == 100. );
	}
	
	/** Issue 61 - Transposing eagerly copies */
	void testIssue61() {
		auto mat1 = Matrix!double([ [1.,2.], [3., 4.] ]);
		auto mat2 = eval( mat1.t );
		Matrix!double mat3;
		mat3[] = mat2.t;
		
		enforce( mat1.cdata == mat2.cdata );
		enforce( mat2.cdata == mat3.cdata );
	}
	
}
