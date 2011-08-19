/** Implementation of the CowArray container, a copy-on-write array which is the default container for vector storage
    types and packed matrix storage types.

    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/
module scid.storage.cowarray;

import scid.internal.assertmessages;
import scid.storage.arraydata;
import scid.blas;
import scid.common.meta, scid.common.storagetraits;
import scid.storage.cowmatrix, scid.matrix;
import scid.ops.expression;
import std.typecons, std.algorithm, std.array;
import std.conv;

version(unittest) {
	import std.stdio;
}

/** A copy-on-write array. Used as a container for storage types. */
struct CowArray( ElementType_ ) {
	alias ElementType_ ElementType;
	alias ArrayData!ElementType Data;
	
	/** Allocate a new array of a given length. Initialize with zero. */
	this()( size_t newLength ) {
		data_.reset( newLength );
		ptr_    = data_.ptr;
		length_ = newLength;
		blas.scal( newLength, Zero!ElementType, ptr_, 1 );
	}
	
	/** Allocate a new uninitialized array of given length. */
	this()( size_t newLength, void* ) {
		data_.reset( newLength );
		ptr_    = data_.ptr;
		length_ = newLength;
	}
	
	/** Allocate a new array and initialize it with the given initializer. */
	this( E )( E[] initializer ) {
		data_.reset( to!(ElementType[])(initializer) );
		ptr_    = data_.ptr;
		length_ = initializer.length;
	}
	
	/** Create an array which is a copy of an existing one. */
	this()( CowArray* other ) {
		data_   = other.data_;
		ptr_    = other.ptr_;
		length_ = other.length_;
	}
	
	/** Resize the array and set all the elements to zero. */
	void resize( size_t len ) {
		resize( len, null );
		blas.scal( len, Zero!ElementType, this.ptr_, 1 );
	}
	
	/** Resize the array and leave the elements uninitialized. */
	void resize( size_t len, void* ) {
		if( len != length || data_.refCount() > 1 ) {
			data_.reset( len );
			ptr_ = data_.ptr;
			length_ = len;
		}
	}
	
	/** Assignment has copy semantics. The actual copy is only performed on modification of the copy however. */
	ref typeof(this) opAssign( CowArray rhs ) {
		swap( rhs.data_, data_ );
		ptr_    = rhs.ptr_;
		length_ = rhs.length_;
		return this;
	}
	
	/** Element access. */
	ElementType index( size_t i ) const
	in {
		checkBounds_( i );
	} body {
		return ptr_[ i ];
	}
	
	/// ditto
	void indexAssign( string op = "" )( ElementType rhs, size_t i )
	in {
		checkBounds_( i );
	} out {
		assert( index( i ) == rhs );
		assert( data_.refCount() == 1 );
	} body {
		unshareData_();
		mixin( "ptr_[ i ]" ~ op ~ "= rhs;" );
	} 
	
	/** Get a slice of this array. */
	typeof( this ) slice( size_t start, size_t end )
	in {
		checkSliceIndices_( start, end );
	} body {
		CowArray r;
		r.data_ = data_;
		r.length_ = end - start;
		r.ptr_    = ptr_ + start;
		return r;
	}
	
	/** Remove the first element. Part of the BidirectionalRange concept. */
	void popFront()
	in {
		checkNotEmpty_!"popFront"();
	} body {
		++ ptr_;
		-- length_;
	}
	
	/** Remove the last element. Part of the BidirectionalRange concept. */
	void popBack()
	in {
		checkNotEmpty_!"popBack"();
	} body {
		-- length_;
	}
	
	@property {
		/** Get a const ref to the wrapped ArrayData. */
		ref const(Data) arrayData() const { return data_; }
		
		/** Get a mutable pointer to the memory used by this storage. */
		ElementType* data() {
			unshareData_();
			return ptr_;
		}
		
		/** Get a const pointer to the memory used by this storage. */
		const(ElementType)* cdata() const {
			return ptr_;
		}
		
		/** Returh the length of this array. Part of the BidirectionalRange concept. */
		size_t length() const {
			return length_;
		}
		
		/** Is the array empty? Part of the BidirectionalRange concept. */
		bool empty() const {
			return length_ == 0;
		}
		
		/** Get the address of this. Needed for a hack to avoid copying in certain cases. */
		typeof(this)* ptr() {
			return &this;
		}
		
		/** Change the first element. Part of the BidirectionalRange concept. */
		void front( ElementType newValue )
		in {
			checkNotEmpty_!"front setter"();
		} body {
			unshareData_();
			*ptr_ = newValue;
		}
		
		/** Change the last element. Part of the BidirectionalRange concept. */
		void back( ElementType newValue  )
		in {
			checkNotEmpty_!"back setter"();
		} body {
			unshareData_();
			*(ptr_ + length_ - 1) = newValue;
		}
		
		/** Return the first element. Part of the BidirectionalRange concept. */
		ElementType front() const
		in {
			checkNotEmpty_!"front"();
		} body {
			return *ptr_;
		}
		
		/** Return the last element. Part of the BidirectionalRange concept. */
		ElementType back() const
		in {
			checkNotEmpty_!"back"();
		} body {
			return *(ptr_ + length_ - 1);
		}
	}
	
	template Promote( T ) {
		static if( isArrayContainer!T || isMatrixContainer!T ) {
			alias CowArrayRef!(Promotion!(BaseElementType!T,ElementType)) Promote;
		} else static if( isScalar!T ) {
			alias CowArrayRef!(Promotion!(T,ElementType)) Promote;
		}
	}
	
private:
	mixin ArrayChecks;
	
	void unshareData_() {
		if( data_.refCount() == 1 )
			return;
		
		// if we only own a slice of data_, we'll only duplicate that
		// bit of the array.
		if( ptr_ == data_.ptr && length_ == data_.length )
			data_.unshare();
		else
			data_.reset( ptr_[ 0 .. length_ ] );
		
		ptr_ = data_.ptr;
	}
	
	Data         data_;
	ElementType* ptr_;
	size_t       length_;
}

/** A simple alias for the preferred reference type for this container. */
template CowArrayRef( T ) {
	alias RefCounted!(CowArray!T, RefCountedAutoInitialize.yes )
			CowArrayRef;
}

unittest {
	// TODO: CowArray tests.
}