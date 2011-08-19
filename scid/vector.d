module scid.vector;

import scid.storage.array;
import scid.storage.arrayview;
import scid.storage.external;
import scid.storage.cowarray;
import scid.storage.constant;
import scid.common.traits;
import scid.common.meta;
import scid.ops.eval;

import std.traits, std.range, std.algorithm, std.conv;

import scid.internal.assertmessages;

enum VectorType {
	Row, Column
}

/**
This template allows for the creation of SciD's vector storage type.  A
SciD Vector is a one-dimensional array-like object that has reference-counted 
copy-on-write semantics.  It can be used in SciD expressions.  VectorType
controls whether the vector is treated as a row or column vector.

Examples:
---
// Create a Vector from an existing array.
double[] arr = [1.0, 2.0, 3.0];
auto vec1 = Vector!double(arr);
vec1[0] = 42;
assert(arr[0] == 1);   // The array is copied, not aliased.
auto vec2 = vec1;
vec2[1] = 84;  
assert(vec1[1] == 2);  // Value semantics.

// Create a Vector from a type and a length.
auto vec3 = Vector!(double, VectorType.Row)(3);
vec3[0] = 8;
vec3[1] = 6;
vec3[2] = 7;

// This is allowed because vec3 is a row vector and vec1 is a column
// vector.
auto expr = vec3 * vec1;

// This is not allowed because vec1 and vec2 are both column vectors.
auto fails = vec1 * vec2;
---
*/
template Vector( ElementOrStorage, VectorType vectorType = VectorType.Column )
		if( isScalar!(BaseElementType!ElementOrStorage) ) {
	
	static if( isScalar!ElementOrStorage )
		alias BasicVector!( ArrayStorage!( ElementOrStorage, vectorType ) ) Vector;
	else
		alias BasicVector!( ElementOrStorage )              Vector;
}

/**
Convenience function for creating a Vector from an existing array, using the
array's element type.

Examples:
---
double[] arr = [1.0, 2.0, 3.0];

// The following are equivalent:
auto vec1 = vector(arr);
auto vec2 = Vector!double(arr);
---
*/
Vector!( T, vectorType ) 
vector( T, VectorType vectorType = VectorType.Column )( T[] array ) {
    return typeof(return)(array);
}

unittest {
    // Test examples, except use convenience functions where possible.
    // Create a Vector from an existing array.
    double[] arr = [1.0, 2.0, 3.0];
    auto vec1 = vector(arr);
    vec1[0] = 42;
    assert(arr[0] == 1);   // The array is copied, not aliased.
    auto vec2 = vec1;
    vec2[1] = 84;  
    assert(vec1[1] == 2);  // Value semantics.

    // Create a Vector from a type and a length.
    auto vec3 = Vector!(double, VectorType.Row)(3);
    vec3[0] = 8;
    vec3[1] = 6;
    vec3[2] = 7;

    // This is allowed because vec3 is a row vector and vec1 is a column
    // vector.
    auto expr = vec3 * vec1;
}

template VectorView( ElementOrStorage, VectorType vectorType = VectorType.Column )
		if( isScalar!(BaseElementType!ElementOrStorage) ) {
			
	alias BasicVector!( ArrayViewStorage!( ElementOrStorage, vectorType ) ) VectorView;
}

template StridedVectorView( ElementOrStorage, VectorType vectorType = VectorType.Column )
		if( isScalar!(BaseElementType!ElementOrStorage) ) {
	
	alias BasicVector!( StridedArrayViewStorage!( ElementOrStorage, vectorType ) ).View StridedVectorView;
}

/**
Template for creating a Vector-like view of a D array.  These have reference
semantics, since they are views instead of full-blown copy-on-write containers.
vectorType controls whether the view is treated by SciD as a row or column
vector.  This can be used in multiple ways:

---
// Create a Vector-like view of an existing array.  This object will have
// reference semantics and will use the storage from the existing array.
double[] arr = [1.0, 2.0, 3.0];
auto view1 = ExternalVectorView!double(arr);
view1[0] = 42;
assert(arr[0] == 42);  // Same storage.
auto view2 = view1;
view2[1] = 84;
assert(view1[1] == 84);  // Reference semantics.
assert(arr[1] == 84);    // Still using the same storage.
---

---
// Create a copy of an existing array using a custom allocator, and obtain
// a vector view of the copy.
auto alloc = newRegionAllocator();
double[] arr = [1.0, 2.0, 3.0];
auto view1 = ExternalVectorView!double(arr, alloc);
view1[0] = 42;
assert(arr[0] == 1);  // Different storage.
auto view2 = view1;
view2[1] = 84;
assert(view1[1] == 84);  // Reference semantics.
assert(arr[1] == 2);    
---

---
// Create a new vector view of length 3 using a custom allocator.
auto alloc = newRegionAllocator();
auto vec1 = ExternalVectorView!double(3, alloc);
vec1[0] = 1;
auto vec2 = vec1;
vec2[0] = 84;
assert(vec1[0] == 84);  // Reference semantics.
---
*/
template ExternalVectorView( ElementOrContainer, VectorType vectorType = VectorType.Column )
		if( isScalar!( BaseElementType!ElementOrContainer ) ) {
	
	static if( isScalar!ElementOrContainer ) {
		alias BasicVector!(
			ArrayViewStorage!(
				ExternalArray!( ElementOrContainer, CowArrayRef!ElementOrContainer ),
				vectorType
			)
		) ExternalVectorView;
	} else {
		alias BasicVector!(
			ArrayViewStorage!(
				ExternalArray!( BaseElementType!ElementOrContainer, ElementOrContainer ),
				vectorType
			)
		) ExternalVectorView;
	}
}

/**
Convenience functions for creating an ExternalVectorView from an array,
using the inferred element type of the array.

Examples:
---
double[] arr = [1.0, 2, 3];

// These two lines are equivalent:
auto view1 = externalVectorView(arr);
auto view2 = ExternalVectorView!(double)(arr);

auto alloc = newRegionAllocator();

// These two lines are also equivalent.
auto view3 = externalVectorView(arr, alloc);
auto view4 = ExternalVectorView!double(arr, alloc);
---
*/
ExternalVectorView!( T, vectorType ) 
externalVectorView( VectorType vectorType = VectorType.Column, T )( T[] array ) {
    return typeof(return)(array);
}

/// Ditto
ExternalVectorView!( T, vectorType ) 
externalVectorView( VectorType vectorType = VectorType.Column, T, Allocator )
( T[] array, Allocator alloc ) {
    return typeof(return)(array, alloc);
}

unittest {
    // Test the examples, except use the convenience functions where possible.
    double[] arr = [1.0, 2.0, 3.0];
    auto view1 = externalVectorView(arr);
    view1[0] = 42;
    assert(arr[0] == 42);  // Same storage.
    auto view2 = view1;
    view2[1] = 84;
    assert(view1[1] == 84);  // Reference semantics.
    assert(arr[1] == 84);    // Still using the same storage.
}

unittest {
    import scid.internal.regionallocator;
    auto alloc = newRegionAllocator();
    double[] arr = [1.0, 2.0, 3.0];
    auto view1 = externalVectorView(arr, alloc);
    view1[0] = 42;
    assert(arr[0] == 1);  // Different storage.
    auto view2 = view1;
    view2[1] = 84;
    assert(view1[1] == 84);  // Reference semantics.
    assert(arr[1] == 2);    
}

unittest {
    import scid.internal.regionallocator;
    auto alloc = newRegionAllocator();
    auto vec1 = ExternalVectorView!double(3, alloc);
    vec1[0] = 1;
    auto vec2 = vec1;
    vec2[0] = 84;
    assert(vec1[0] == 84);  // Reference semantics.
}

struct BasicVector( Storage_ ) {
	alias BaseElementType!Storage                          ElementType;
	alias Storage_                                         Storage;
	alias BasicVector!( typeof(Storage.init.view(0,0)) )   View;
	alias BasicVector!( Storage.Transposed )               Transposed;
	alias storage                                          this;
	
	static if( is( Storage.Temporary ) )
		alias BasicVector!( Storage.Temporary ) Temporary;
	else
		alias typeof( this ) Temporary;
	
	static if( is( typeof(Storage.vectorType) ) )
		alias Storage.vectorType vectorType;
	else
		alias VectorType.Column vectorType;
	
	/** Whether the storage can be resized. */
	enum isResizable = is( typeof( Storage.init.resize(0) ) );
	
	static if( is( typeof(Storage.init.view(0,0,0)) R ) )
		alias BasicVector!R StridedView;
	
	//static assert( isVectorStorage!Storage );

	this( A... )( A args ) if( A.length > 0 && !is( A[0] : Storage ) && !isVector!(A[0]) && !isExpression!(A[0]) ) {
		storage = Storage(args);
	}
	
	this( Expr )( Expr expr ) if( isExpression!Expr ) {
		this[] = expr;
	}
	
	this( A )( BasicVector!A other ) {
		static if( is( A : Storage ) ) storage = other.storage;
		else                           this[] = other;
	}
	
	this()( auto ref Storage stor ) {
		storage = stor;
	}
	
	ElementType opIndex( size_t i ) const {
		return storage.index( i );
	}
	
	void opIndexAssign( ElementType rhs, size_t i ) {
		storage.indexAssign( rhs, i );
	}
	
	void opIndexOpAssign( string op )( ElementType rhs, size_t i ) {
		storage.indexAssign!op( rhs, i );
	}
	
	ref typeof(this) opAssign( typeof(this) rhs ) {
		move( rhs.storage, storage );
		return this;
	}
	
	typeof( this ) opSlice() {
		return typeof(this)( storage );
	}
	
	typeof( this ) opSlice( size_t start, size_t end ) {
		return typeof(this)( storage.slice( start, end ) );	
	}
	
	bool opEquals( Rhs )( Rhs rhs ) const if( isInputRange!Rhs ) {
		size_t i = 0;
		foreach( x ; rhs ) {
			if( i >= length || this[ i ++ ] != x )
				return false;
		}
		
		return true;
	}
	
	/** Resize the vector and leave the memory uninitialized. If not resizeable simply check that the length is
	    correct.
	*/
	void resize( size_t newLength, void* ) {
		static if( isResizable ) {
			storage.resize( newLength, null );
		} else {
			assert( length == newLength,
				lengthMismatch_( newLength ) );
		}
	}
	
	/** Resize the vectors and set all the elements to zero. If not resizeable, check that the length is correct
	    and just set the elements to zero.
	*/
	void resize( size_t newLength ) {
		static if( isResizable ) {
			storage.resize( newLength );
		} else {
			this.resize( newLength, null );
			evalScaling( Zero!ElementType, this );
		}
	}

	void opSliceAssign( Rhs )( auto ref Rhs rhs ) {
		static if( is( Rhs E : E[] ) && isConvertible!( E, ElementType  ) )
			evalCopy( ExternalVectorView!( E, vectorType )( rhs ), this );
		else static if( closureOf!Rhs == Closure.Scalar )
			evalCopy( relatedConstant( rhs, this ), this );
		else {
			evalCopy( rhs, this );
		}
		
	}
	
	void opSliceAssign( Rhs )( Rhs rhs, size_t start, size_t end ) {
		view( start, end )[] = rhs;
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs ) if( op == "+" || op == "-" ) {
		enum scalarRhs = closureOf!Rhs == Closure.Scalar;
		static if( op == "+" ) {
			static if( scalarRhs ) evalScaledAddition( One!ElementType, relatedConstant(rhs, this), this );
			else                   evalScaledAddition( One!ElementType, rhs, this );
		} else static if( op == "-" ) {
			static if( scalarRhs ) evalScaledAddition( One!ElementType, relatedConstant(-rhs, this), this );
			else                   evalScaledAddition( MinusOne!ElementType, rhs, this );
		}
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs ) if( (op == "*" || op == "/") && isConvertible!(Rhs,ElementType) ) {
		static if( op == "/" )
			rhs = One!ElementType / rhs;
		evalScaling( to!ElementType(rhs), this );
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs, size_t start, size_t end ) {
		mixin( "view( start, end )[] " ~ op ~ "= rhs;" );
	}
	
	View view( size_t start, size_t end ) {
		return typeof( return )( storage.view( start, end ) );
	}
	
	static if( is( StridedView ) ) {
		StridedView view( size_t start, size_t end, size_t stride ) {
			return typeof( return )( storage.view( start, end, stride ) );	
		}
	}
	
	static if( isInputRange!Storage ) {
		void popFront() { storage.popFront(); }
		void popBack()  { storage.popBack(); }
	}
	
	@property {
		bool empty() const {
			static if( is(typeof(Storage.init.empty)) )
				return storage.empty;
			else
				return storage.length == 0;
		}
		
		size_t length() const {
			return storage.length;
		}
		
		ElementType front() const
		in {
			checkNotEmpty_!"front"();
		} body {
			static if( is( typeof(Storage.init.front) ) )
				return storage.front;
			else
				return storage.index( 0 );
		}
		
		ElementType back() const
		in {
			checkNotEmpty_!"front"();
		} body {
			static if( is( typeof(Storage.init.back) ) )
				return storage.back;
			else
				return storage.index( storage.length - 1 );
		}
	}
	
	string toString() const {
		if( empty )
			return "[]";
		
		auto r = appender!string("[");
		r.put( to!string( this[ 0 ] ) );
		foreach( i ; 1 .. length ) {
			r.put( ", " );
			r.put( to!string( this[i] ) );
		}
		r.put( "]" );
		return r.data();
	}
	
	alias toString pretty;
	
	static if( vectorType == VectorType.Column ) mixin Operand!( Closure.ColumnVector );
	else                                         mixin Operand!( Closure.RowVector    );
	
	template Promote( T ) {
		static if( isVector!T ) {
			alias BasicVector!( Promotion!(Storage,T.Storage) ) Promote;
		} else static if( isMatrix!T ) {
			alias BasicVector!( Promotion!(Storage,T.Storage) ) Promote;
		} else static if( isScalar!T ) {
			alias BasicVector!( Promotion!(Storage,T) ) Promote;
		}
	}
	
	Storage storage;
	
private:
	mixin ArrayChecks;
}

unittest {
	// TODO: Write tests for Vector.
}
