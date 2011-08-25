/** Common error messages to use with assertions.

    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/
module scid.internal.assertmessages;

import std.string;
import std.range, std.conv;

private alias to!string tos;

/** Checks and error messages for array-like structs. */
mixin template ArrayChecks() {
	/** Check that an index is within the bounds. */
	void checkBounds_( size_t i ) const {
		auto len = this.length;
		assert( i < len, "Out of bounds: i = " ~ tos( i ) ~ ", length = " ~ tos( len ) );
	}
	
	/** Check that slice indices are valid. */
	void checkSliceIndices_( size_t start, size_t end ) const {
		auto len = this.length;
		assert( start <= end && end <= len,
			"Invalid slicing indices: start = " ~ tos( start ) ~ ", end = " ~ tos( end ) );
	}
	
	/** Check that the lengths in a slice assignment match. */
	void checkSliceAssignLength_( size_t start, size_t end, size_t rhsLength ) const {
		auto sliceLength = end - start;
		auto len = this.length;
		assert( start <= end && end <= len && sliceLength == rhsLength,
			"Length mismatch in slice assignment: start = " ~ tos( start ) ~
			", end = " ~ tos( end ) ~
			", rhsLength = " ~ tos( rhsLength )
		);
	}
	
	/** Check that the lengths in an assignment match. */
	void checkAssignLength_( size_t rhsLength ) const {
		auto len = this.length;
		assert( len == rhsLength,
			"Length mismatch in assignment: length = " ~ tos( len ) ~ ", rhsLength = " ~ tos( rhsLength ) ~ "."
		);
	}
	
	/** Check that the range is not empty. */
	void checkNotEmpty_( string op = "function" )() const {
		assert( !this.empty,
			"Invalid " ~ op ~ " call on empty range." );
	}
}

/** Checks and error messages for matrix-like structs (matrix literal, storages and containers). */
mixin template MatrixChecks() {
	/** Check that an index is within the bounds. */
	void checkBounds_( size_t i, size_t j ) const {
		auto r = this.rows, c = this.columns;
		assert( i < r && j < c, "Out of bounds: (i,j) = (" ~ tos( i ) ~ ", " ~ tos( j ) ~ "), dimensions = (" ~ tos( r ) ~ ", " ~ tos( c ) ~ ")" );
	}
	
	/** Check that slice indices are valid. */
	void checkSliceIndices_( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) const {
		auto r = this.rows, c = this.columns;
		assert( rowStart <= rowEnd && rowEnd <= r && colStart <= colEnd && colEnd <= c,
			"Invalid slicing indices [ " ~ tos( rowStart ) ~ " .. " ~ tos( rowEnd ) ~ " ][ "
			~ tos( colStart ) ~ " .. " ~ tos( colEnd ) ~ "], dimensions = (" ~ tos( r ) ~ ", " ~ tos( c ) ~ ")"
		);
	}
	
	/** Check that a built-in array of arrays is a valid general matrix. */
	static void checkGeneralInitializer_( E )( E[][] mat ) {
		debug {
			if( !mat.length )
				return;
			
			size_t cols = mat[ 0 ].length;
			foreach( i, row ; mat )
				assert( row.length == cols,
					"Inconsistent number of columns in general matrix initializer: row(0).length = " ~ tos( cols ) ~
					", row(" ~ tos( i ) ~ ").length = " ~ tos( row.length )
				);
		}
	}
	
	/** Check that a major dimension and an array form a valid general matrix initializer. */
	static void checkGeneralInitializer_( E )( size_t majorDimension, E[] initializer ) {
		assert( ( (majorDimension != 0) ^ initializer.empty) && initializer.length % majorDimension == 0,
			"Invalid initializer for general matrix: majorDimension = " ~ tos( majorDimension ) ~
			", initializer.length = " ~ tos( initializer.length )
		);
	}
	
	/** Check triangular initializer. */
	static void checkTriangularInitializer_( E )( E[] initializer ) {
		debug {
			import std.math;
			
			auto tri  = (sqrt( initializer.length * 8.0 + 1.0 ) - 1.0 ) / 2.0;
			assert( tri - cast(int) tri <= 0,
				"Initializer length is not triangular number: initializer.length = " ~ tos( initializer.length )
			);
		}
	}
	
	/** Check that dimensions passed to a square matrix type are equal. */
	static void checkSquareDims_( string matrixType = "square" )( size_t newRows, size_t newCols ) {
		assert( newRows == newCols,
				"Non-square dimensions for " ~ matrixType ~ " matrix (" ~ tos( newRows ) ~ ", " ~ tos( newCols ) ~ ")"
		);
	}
	
	/** Check that the range is not empty. */
	void checkNotEmpty_( string op = "function" )() const {
		assert( !this.empty,
			"Invalid " ~ op ~ " call on empty range." );
	}
	
	/** Check that dimensions match in assignment. */
	void checkAssignDims_( size_t rhsRows, size_t rhsColumns ) const {
		auto r = this.rows, c = this.columns;
		assert( r == rhsRows && c == rhsColumns,
			"Dimension mismatch in matrix assignment: lhsDims = (" ~ tos( r ) ~ ", " ~ tos( c ) ~ "), rhsDims = (" ~
			tos( rhsRows ) ~ ", " ~ tos( rhsColumns ) ~ ")"
		);
	}
}

string stridedToString( S )( const(S)* ptr, size_t len, size_t stride ) {
		if( len == 0 )
			return "[]";
	
		auto app = appender!string("[");
		app.put( to!string(*ptr) );
		auto e = ptr + len * stride;
		for( ptr += stride; ptr < e ; ptr += stride ) {
			app.put( ", " );
			app.put( to!string( *ptr ) );
		} 
		app.put(']');
		
		return app.data();
	}
	
string matrixToString( S )( char trans, size_t m, size_t n, const(S)* a, size_t lda )
in {
	assert( a || (!m || !n) );
} body {
	if( m == 0 || n == 0 )
		return "[]";
		
	auto app = appender!string("[");
	if( trans == 'n' ) {
		app.put( stridedToString( a, n, lda ) );
		auto e = a + m;
		for( ++ a; a < e ; ++a ) {
			app.put( ", " );
			app.put( stridedToString( a, n, lda ) );
		}
		app.put(']');
	} else {
		app.put( stridedToString(a, m, 1) );
		auto e = a + n * lda;
		for( a += lda ; a < e ; a += lda ) {
			app.put( ", " );
			app.put( stridedToString( a, m, 1 ) );
		}
		app.put(']');
	}
		
	return app.data();
}