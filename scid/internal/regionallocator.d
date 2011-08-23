/**
$(D RegionAllocator) is a memory allocator based on segmented
stacks.  A segmented stack is similar to a regular stack in that memory is
allocated and freed in last in, first out order.  When memory is requested from
a segmented stack, it first checks whether enough space is available in the
current segment, and if so increments the stack pointer and returns.  If not,
a new segment is allocated.  When memory is freed, the stack pointer is
decremented.  If the last segment is no longer in use, it may be returned to
where it was allocated from or retained for future use.

$(D RegionAllocator) has the following advantages compared to allocation on the
call stack:

1.  Pointers to memory allocated on a $(D RegionAllocator) stack are still
    valid when the function they were allocated from returns, unless the
    last instance of the RegionAllocator object they were allocated from
    goes out of scope.  Functions can be written to create and return data
    structures on the $(D RegionAllocator) stack.
    
2.  It's possible to have more than one stack per thread.

3.  Since it is a segmented stack, large allocations can be performed with no
    danger of stack overflow errors.

4.  It guarantees 16-byte alignment of all allocated objects.

It has the following disadvantages:

1.  The extra checks due to stack segmentation make allocation slightly
    slower and prevent simultanelous deallocation of several objects from
    being a single decrement of the stack pointer.

2.  Memory allocated on RegionAllocator is accessed through an extra layer of
    pointer indirection compared to memory on the call stack.

It has the following advantages compared to heap allocation:

1.  Both allocation and deallocation are extremely fast.  Most allocations
    consist of verifying enough space is available, incrementing a pointer and
    a performing a few cheap bookkeeping operations.  Most deallocations
    consist decrementing a pointer and performing a few cheap bookkeeping
    operations.

2.  The segmented stack is thread-local, so synchronization is only needed
    when a segment needs to be allocated or freed.

3.  Fragmentation is not an issue when allocating memory on the
    $(D RegionAllocator) stack, though it can be an issue when trying to allocate
    a new segment.

It has the following disadvantages compared to heap allocation:

1.  The requirement that memory be freed in last in, first out order.

2.  No automatic garbage collection.

A segmented stack may be created manually.  Alternatively, a default 
thread-local stack that is automatically created lazily on the first
attempt to use it may be used.

Synopsis:
---
void fun() {
    // Create a new RegionAllocator using the default thread-local stack.
    auto alloc = newRegionAllocator();

    // Allocate a temporary array on the RegionAllocator stack.
    auto arr = alloc.uninitializedArray!(double[][][])(8, 6, 7);
    assert(arr.length == 8);
    assert(arr[0].length == 6);
    assert(arr[0][0].length == 7);

    // When alloc goes out of scope, all memory allocated by it is freed.
    // If alloc has been copied, then the memory is freed when all copies
    // have gone out of scope.
}
---

Author:  David Simcha
Copyright:  Copyright (c) 2008-2011, David Simcha.
License:    $(WEB boost.org/LICENSE_1_0.txt, Boost License 1.0)
*/

module scid.internal.regionallocator;

import std.traits, core.memory, std.range, core.exception, std.conv,
    std.algorithm, std.typetuple, std.exception, std.typecons;

static import core.stdc.stdlib;

// This is just for convenience/code readability/saving typing.
private enum ptrSize = (void*).sizeof;

// This was accidentally assumed in a few places and I'm too lazy to fix it
// until I see proof that it needs to be fixed.
static assert(bool.sizeof == 1);

enum size_t defaultSegmentSize = 4 * 1_024 * 1_024;

/**
The exception that is thrown on invalid use of $(RegionAllocator) and
$(D RegionAllocatorStack).  This exception is not thrown on out of memory.
An $(D OutOfMemoryError) is thrown instead.
*/
class RegionAllocatorException : Exception {
    this(string msg, string file, int line) {
        super(msg, file, line);
    }
}

/**
This flag determines whether a given $(D RegionAllocatorStack) is scanned for 
pointers by the garbage collector (GC).  If yes, the entire stack is scanned, 
not just the part currently in use, since there is currently no efficient way to 
modify the bounds of a GC region.  The stack is scanned conservatively, meaning 
that any bit pattern that would point to GC-allocated memory if interpreted as 
a pointer is considered to be a pointer.  This can result in GC-allocated
memory being retained when it should be freed.  Due to these caveats,
it is recommended that any stack scanned by the GC be small and/or short-lived.
*/
enum GCScan : bool {
    ///
    no = false,
    
    ///
    yes = true
}

/**
This object represents a segmented stack.  Memory can be allocated from this
stack using a $(XREF regionallocator RegionAllocator) object.  Multiple 
$(D RegionAllocator) objects may be created per 
$(D RegionAllocatorStack) but each $(D RegionAllocator) uses a single 
$(D RegionAllocatorStack). 

For most use cases it's convenient to use the default thread-local
instance of $(D RegionAllocatorStack), which is lazily instantiated on
the first call to the global function 
$(XREF regionallocator, newRegionAllocator).  Occasionally it may be useful
to have multiple independent stacks in one thread, in which case a 
$(D RegionAllocatorStack) can be created manually.

$(D RegionAllocatorStack) is reference counted and has reference semantics.
When the last copy of a given instance goes out of scope, the memory 
held by the $(D RegionAllocatorStack) instance is released back to the
heap.  This cannot happen before memory allocated to a $(D RegionAllocator)
instance is released back to the stack, because a $(D RegionAllocator)
holds a copy of the $(D RegionAllocatorStack) instance it uses.

Examples:
---
import std.regionallocator;

void main() {
    fun1();
}

void fun1() {
    auto stack = RegionAllocatorStack(1_048_576, GCScan.no);
    fun2(stack);
    
    // At the end of fun1, the last copy of the RegionAllocatorStack
    // instance pointed to by stack goes out of scope.  The memory
    // held by stack is released back to the heap.
}

void fun2(RegionAllocatorStack stack) {
    auto alloc = stack.newRegionAllocator();
    auto arr = alloc.newArray!(double[])(1_024);
    
    // At the end of fun2, the last copy of the RegionAllocator instance
    // pointed to by alloc goes out of scope.  The memory used by arr 
    // is released back to stack.
}
---
*/
struct RegionAllocatorStack {
private:
    RefCounted!(RegionAllocatorStackImpl, RefCountedAutoInitialize.no) impl;
    bool initialized;
    bool _gcScanned;

public:    
    /**
    Create a new $(D RegionAllocatorStack) with a given segment size in bytes.
    */
    this(size_t segmentSize, GCScan shouldScan) {
        this._gcScanned = shouldScan;
        if(segmentSize == 0) {
            throw new RegionAllocatorException( 
                "Cannot create a RegionAllocatorStack with segment size of 0.",
                __FILE__, __LINE__
            );
        }
        
        impl = typeof(impl)(segmentSize, shouldScan);
        initialized = true;
    }
    
    /**
    Creates a new $(D RegionAllocator) region using this stack.  
    */    
    RegionAllocator newRegionAllocator() {
        if(!initialized) {
            throw new RegionAllocatorException(
                "Cannot create a RegionAllocator from an " ~
                "uninitialized RegionAllocatorStack.  Did you call " ~
                "RegionAllocatorStack's constructor?",
                __FILE__,
                __LINE__
            );
        }
        
        return RegionAllocator(this);
    }
    
    /**
    Whether this stack is scanned by the garbage collector.
    */
    bool gcScanned() @property const pure nothrow @safe {
        return _gcScanned;
    }
}

private struct RegionAllocatorStackImpl {
   
    this(size_t segmentSize, GCScan shouldScan) {
        this.segmentSize = segmentSize;
        space = alignedMalloc(segmentSize, shouldScan);

        // We don't need 16-byte alignment for the bookkeeping array.
        immutable nBookKeep = segmentSize / alignBytes;
        bookkeep = (cast(SizetPtr*) core.stdc.stdlib.malloc(nBookKeep))
                    [0..nBookKeep / SizetPtr.sizeof];
        
        if(!bookkeep.ptr) {
            outOfMemory();
        }
        
        nBlocks++;
    }
    
    size_t segmentSize;  // The size of each segment.
    
    size_t used;
    void* space;
    size_t bookkeepIndex;
    SizetPtr[] bookkeep;
    uint nBlocks;
    uint nFree;
    size_t regionIndex = size_t.max;

    // inUse holds info for all blocks except the one currently being
    // allocated from.  freeList holds space ptrs for all free blocks.
    
    static struct Block {
        size_t used = 0;
        void* space = null;
    }
    
    SimpleStack!(Block) inUse;
    SimpleStack!(void*) freeList;

    void doubleSize(ref SizetPtr[] bookkeep) {
        size_t newSize = bookkeep.length * 2;
        auto ptr = cast(SizetPtr*) core.stdc.stdlib.realloc(
            bookkeep.ptr, newSize * SizetPtr.sizeof);
            
        if(!ptr) {
            outOfMemory();
        }
        
        bookkeep = ptr[0..newSize];
    }

    // Add an element to bookkeep, checking length first.
    void putLast(void* last) {
        if (bookkeepIndex == bookkeep.length) doubleSize(bookkeep);
        bookkeep[bookkeepIndex].p = last;
        bookkeepIndex++;
    }

    void putLast(size_t num) {
        return putLast(cast(void*) num);
    }
    
    // The number of objects allocated on the C heap instead of here.
    ref size_t nLargeObjects() @property pure nothrow {
        return bookkeep[regionIndex - 4].i;
    }
    
    // The number of extra segments that have been allocated in the current
    // region beyond the base segment of the region.
    ref size_t nExtraSegments() @property pure nothrow {
        return bookkeep[regionIndex - 3].i;
    }
    
    // Pushes a segment to the internal free list and frees segments to the
    // heap if there are more than 2x as many segments on the free list as
    // in use.
    void freeSegment() {
        nExtraSegments--;
        freeList.push(space);
        auto newHead = inUse.pop();
        space = newHead.space;
        used = newHead.used;
        nBlocks--;
        nFree++;

        if (nFree >= nBlocks * 2) {
            foreach(i; 0..nFree / 2) {
                alignedFree(freeList.pop);
                nFree--;
            }
        }
    }
    
    void destroy() {
        if(space) {
            alignedFree(space);
            space = null;
        }

        if(bookkeep) {
            core.stdc.stdlib.free(bookkeep.ptr);
            bookkeep = null;
        }

        while(inUse.index > 0) {
            auto toFree = inUse.pop();
            alignedFree(toFree.space);
        }

        inUse.destroy();

        while(freeList.index > 0) {
            auto toFree = freeList.pop();
            alignedFree(toFree);
        }

        freeList.destroy();
    }

    ~this() {
        destroy();
    }
}

/**
These properties get and set the segment size of the default thread-local
$(D RegionAllocatorStack) instance.  The default size is 4 megabytes.
The setter is only effective before the global function
$(D newRegionAllocator) has been called for the first time in the current
thread.  Attempts to set this property after the first call to this
function from the current thread throw a $(D RegionAllocatorException).
*/
size_t threadLocalSegmentSize() @property nothrow @safe {
    return _threadLocalSegmentSize;
}

/// Ditto
size_t threadLocalSegmentSize(size_t newSize) @property @safe {
    if(threadLocalInitialized) {
        throw new RegionAllocatorException(
            "Cannot set threadLocalSegmentSize after the thread-local " ~
            "RegionAllocatorStack has been used for the first time.",
            __FILE__,
            __LINE__
        );
    }
    
    return _threadLocalSegmentSize = newSize;
}

/**
These properties determine whether the default thread-local 
$(D RegionAllocatorStack) instance is scanned by the garbage collector.
The default is no.  In most cases, scanning a stack this long-lived is not
recommended, as it will cause too many false pointers.  (See $(XREF 
regionallocator, GCScan) for details.)  

The setter is only effective before the global function
$(D newRegionAllocator) has been called for the first time in the current
thread.  Attempts to set this property after the first call to this
function from the current thread throw a $(D RegionAllocatorException).
*/
bool scanThreadLocalStack() @property nothrow @safe {
    return _scanThreadLocalStack;
}

/// Ditto
bool scanThreadLocalStack(bool shouldScan) @property @safe {
    if(threadLocalInitialized) {
        throw new RegionAllocatorException(
            "Cannot set scanThreadLocalStack after the thread-local " ~
            "RegionAllocatorStack has been used for the first time.",
            __FILE__,
            __LINE__
        );
    }
    
    return _scanThreadLocalStack = shouldScan;
}

private size_t _threadLocalSegmentSize = defaultSegmentSize;
private RegionAllocatorStack threadLocalStack;
private bool threadLocalInitialized;
private bool _scanThreadLocalStack = false;

// Ensures the thread-local stack is initialized, then returns it.
private ref RegionAllocatorStack getThreadLocal() {
    if(!threadLocalInitialized) {
        threadLocalInitialized = true;
        threadLocalStack = RegionAllocatorStack(
            threadLocalSegmentSize, cast(GCScan) scanThreadLocalStack
        );
    }
    
    return threadLocalStack;
}

static ~this() {
    if(threadLocalInitialized) {
        threadLocalStack.impl.refCountedPayload.destroy();
    }
}

/**
This struct provides an interface to the $(D RegionAllocator) functionality
and enforces scoped deletion.  A new instance using the thread-local 
$(D RegionAllocatorStack) instance is created using the global
$(D newRegionAllocator) function.  A new instance using 
an explicitly created $(D RegionAllocatorStack) is created using 
$(D RegionAllocatorStack.newRegionAllocator).

Each instance has reference semantics in that any copy will allocate from the 
same memory.  When the last copy of an instance goes out of scope, all memory 
allocated via that instance is freed.  Only the most recently created 
still-existing $(D RegionAllocator) using a given $(D RegionAllocatorStack)  
may be used to allocate and free memory at any given time.  Deviations
from this model result in a $(D RegionAllocatorException) being thrown.

An uninitialized $(D RegionAllocator) (for example $(D RegionAllocator.init))
has semantics similar to a null pointer.  It may be assigned to or passed to
a function.  However, any attempt to call a method will result in a
$(D RegionAllocatorException) being thrown.  

Examples:
---
void foo() {
    auto alloc = newRegionAllocator();
    auto ptr1 = bar(alloc);
    auto ptr2 = alloc.allocate(42);

    // The last copy of the RegionAllocator object used to allocate ptr1
    // and ptr2 is going out of scope here.  The memory pointed to by
    // both ptr1 and ptr2 will be freed.
}

void* bar(RegionAllocator alloc) {
    auto ret = alloc.allocate(42);

    auto alloc2 = newRegionAllocator();
    auto ptr3 = alloc2.allocate(42);

    // ptr3 was allocated using alloc2, which is going out of scope.
    // Its memory will therefore be freed.  ret was allocated using alloc.
    // A copy of this RegionAllocator is still alive in foo() after
    // bar() executes.  Therefore, ret will not be freed on returning and
    // is still valid after bar() returns.

    return ret;
}

void* thisIsSafe() {
    // This is safe because the two RegionAllocator objects being used
    // are using two different RegionAllocatorStack objects.
    auto alloc = newRegionAllocator();
    auto ptr1 = alloc.allocate(42);
    
    auto stack = RegionAllocatorStack(1_048_576, GCScan.no);
    auto alloc2 = stack.newRegionAllocator();
    
    auto ptr2 = alloc2.allocate(42);
    auto ptr3 = alloc.allocate(42);
}    

void* dontDoThis() {
    auto alloc = newRegionAllocator();
    auto ptr1 = alloc.allocate(42);
    auto alloc2 = newRegionAllocator();

    // Error:  Allocating from a RegionAllocator instance other than the
    // most recently created one that's still alive from a given stack.
    auto ptr = alloc.allocate(42);
}

void uninitialized() {
    RegionAllocator alloc;
    auto ptr = alloc.allocate(42);  // Error:  alloc is not initialized.
    auto alloc2 = alloc;  // Ok.  Both alloc, alloc2 are uninitialized.
    
    alloc2 = newRegionAllocator();
    auto ptr2 = alloc2.allocate(42);  // Ok.
    auto ptr3 = alloc.allocate(42);  // Error:  alloc is still uninitialized.
    
    alloc = alloc2;
    auto ptr4 = alloc.allocate(42);  // Ok.
}    
---

Note:  Allocations larger than $(D this.segmentSize) are handled as a special
case and fall back to allocating directly from the C heap.  These large 
allocations are freed as if they were allocated on a $(D RegionAllocatorStack) 
when $(D free) or $(D freeLast) is called or the last copy of a 
$(D RegionAllocator) instance goes out of scope.  However, due to the extra 
bookkeeping required, destroying a region (as happens when the last copy of
a $(D RegionAllocator) instance goes out of scope) will require time linear
instead of constant in the number of allocations for regions where these
large allocations are present.
*/
struct RegionAllocator {
private:
    RegionAllocatorStack stack;
    
    // The region index that should be current anytime this instance is
    // being used.  This is checked for in allocate() and free().
    size_t correctRegionIndex = size_t.max;
    
    this(ref RegionAllocatorStack stack) {
        assert(stack.initialized);
        auto impl = &(stack.impl.refCountedPayload());
        this.stack = stack;
        
        with(*impl) {
            putLast(0);            // nLargeObjects.
            putLast(0);            // nExtraSegments.
            putLast(regionIndex);  // Old regionIndex.
            putLast(1);            // Ref count of current RegionAllocator.
            regionIndex = bookkeepIndex;
            correctRegionIndex = regionIndex;
        }
    }
    
    // CTFE function, for static assertions.  Can't use bsr/bsf b/c it has
    // to be usable at compile time.
    static bool isPowerOf2(size_t num) pure nothrow @safe {
        return num && !(num & (num - 1));
    }
    
    alias RegionAllocatorStackImpl Impl;  // Save typing.

    // This is written as a mixin instead of a function because it's 
    // performance critical and it wouldn't be inlinable because it throws.
    // By using a mixin, initialized can be checked all the time instead of
    // just in debug mode, for negligible performance cost.
    enum string getStackImplMixin = q{
        if(!initialized) {
            throw new RegionAllocatorException(
                "RegionAllocator instance not initialized.  Please use " ~
                "newRegionAllocator() to create a RegionAllocator object.",
                __FILE__,
                __LINE__
            );
        }
        
        auto impl = &(stack.impl.refCountedPayload());
    };

    void incrementRefCount() {
        mixin(getStackImplMixin);
        impl.bookkeep[correctRegionIndex - 1].i++;
    }

    void decrementRefCount() {
        mixin(getStackImplMixin);
        impl.bookkeep[correctRegionIndex - 1].i--;
        
        if(impl.bookkeep[correctRegionIndex - 1].i == 0) {
            if(impl.regionIndex != correctRegionIndex) {
                throw new RegionAllocatorException(
                    "Cannot free RegionAlloc regions in non-last in first " ~
                    "out order.  Did you return a RegionAllocator from a " ~
                    "function?", 
                    __FILE__,
                    __LINE__
                );
            }

            with(*impl) {
                // Free allocations one at a time until we don't have any
                // more large objects.
                while (nLargeObjects > 0 && bookkeepIndex > regionIndex) {
                    assert(bookkeepIndex > regionIndex);
                    freeLast();
                }
                
                // Free any extra segments that were used by this region.
                while(nExtraSegments > 0) {
                    freeSegment();
                }
                
                if(bookkeepIndex > regionIndex) {
                    // Then there's something left to free.
                    used = bookkeep[regionIndex].i - cast(size_t) space;
                }
                
                bookkeepIndex = regionIndex - 4;
                regionIndex = bookkeep[regionIndex - 2].i;
            }
        }
    }

    bool initialized() @property const pure nothrow @safe {
        return correctRegionIndex < size_t.max;
    }
    
    Unqual!(ElementType!(R))[] arrayImplStack(R)(R range) {
        alias ElementType!(R) E;
        alias Unqual!(E) U;
        static if(hasLength!(R)) {
            U[] ret = uninitializedArray!(U[])(range.length);
            copy(range, ret);
            return ret;
        } else {
            mixin(getStackImplMixin);
            auto startPtr = allocate(0);
            size_t bytesCopied = 0;

            while(!range.empty) {
                auto elem = range.front;
                if(impl.used + U.sizeof <= segmentSize) {
                    range.popFront;
                    *(cast(U*) (startPtr + bytesCopied)) = elem;
                    bytesCopied += U.sizeof;
                    impl.used += U.sizeof;
                } else {
                    if(bytesCopied + U.sizeof >= segmentSize / 2) {
                        // Then just heap-allocate.
                        U[] result = (cast(U*) 
                            alignedMalloc(bytesCopied * 2, gcScanned))
                            [0..bytesCopied / U.sizeof * 2];

                        immutable elemsCopied = bytesCopied / U.sizeof;
                        result[0..elemsCopied] = (cast(U*) startPtr)
                            [0..elemsCopied];
                        finishCopy(result, range, elemsCopied);
                        freeLast();
                        impl.putLast(result.ptr);
                        impl.nLargeObjects++;
                        return result;
                    } else {
                        // Force allocation of a new segment.
                        U[] oldData = (cast(U*) startPtr)
                            [0..bytesCopied / U.sizeof];
                        impl.used -= bytesCopied;
                        impl.bookkeepIndex--;
                        U[] arr = uninitializedArray!(U[])
                            (bytesCopied / U.sizeof + 1);
                        arr[0..oldData.length] = oldData[];
                        startPtr = impl.space;
                        arr[$ - 1] = elem;
                        bytesCopied += U.sizeof;
                        range.popFront;
                    }
                }
            }
            auto rem = bytesCopied % .alignBytes;
            if(rem != 0) {
                auto toAdd = .alignBytes - rem;
                if(impl.used + toAdd < RegionAllocator.segmentSize) {
                    impl.used += toAdd;
                } else {
                    impl.used = RegionAllocator.segmentSize;
                }
            }
            return (cast(U*) startPtr)[0..bytesCopied / U.sizeof];
        }
    }

    Unqual!(ElementType!(R))[] arrayImplHeap(R)(R range) {
        // Initial guess of how much space to allocate.  It's relatively large 
        // b/c the object will be short lived, so speed is more important than 
        // space efficiency.
        enum initialGuess = 128;

        mixin(getStackImplMixin);
        alias Unqual!(ElementType!R) E;
        auto arr = (cast(E*) alignedMalloc(E.sizeof * initialGuess, true))
            [0..initialGuess];

        finishCopy(arr, range, 0);
        impl.putLast(arr.ptr);
        impl.nLargeObjects++;
        return arr;
    }
    
public:
    
    this(this) {
        if(initialized) incrementRefCount();
    }

    ~this() {
        if(initialized) decrementRefCount();
    }

    void opAssign(RegionAllocator rhs) {
        if(initialized) decrementRefCount();
        this.stack = rhs.stack;
        this.correctRegionIndex = rhs.correctRegionIndex;
        if(initialized) incrementRefCount();
    }

    /**
    Allocates $(D nBytes) bytes on the $(D RegionAllocatorStack) used by this
    $(D RegionAllocator) instance.  The last block allocated from this 
    $(D RegionAllocator) instance can be freed by calling
    $(D RegionAllocator.free) or $(D RegionAllocator.freeLast) or will be
    automatically freed when the last copy of this $(D RegionAllocator)
    instance goes out of scope.
    
    Allocation requests larger than $(D segmentSize) are
    allocated directly on the C heap, are scanned by the GC iff
    the $(D RegionAllocatorStack) instance that this object uses is scanned by
    the GC, and are freed according to the same rules as described above.
    */
    void* allocate(size_t nBytes) {
        mixin(getStackImplMixin);
        if(impl.regionIndex != this.correctRegionIndex) {
            throw new RegionAllocatorException(
                "Cannot allocate memory from a RegionAllocator that is not " ~
                "currently at the top of the stack.",
                __FILE__,
                __LINE__
            );
        }
        
        nBytes = allocSize(nBytes);
        with(*impl) {
            void* ret;
            if (segmentSize - used >= nBytes) {
                ret = space + used;
                used += nBytes;
            } else if (nBytes > segmentSize) {
                ret = alignedMalloc(nBytes, gcScanned);
                impl.nLargeObjects++;
            } else if (nFree > 0) {
                nExtraSegments++;
                inUse.push(Block(used, space));
                space = freeList.pop;
                used = nBytes;
                nFree--;
                nBlocks++;
                ret = space;
            } else { // Allocate more space.
                nExtraSegments++;
                inUse.push(Block(used, space));
                space = alignedMalloc(segmentSize, gcScanned);
                nBlocks++;
                used = nBytes;
                ret = space;
            }
            putLast(ret);
            return ret;
        }
    }

    /**
    Frees the last block of memory allocated by the current
    $(D RegionAllocator).  Throws a $(D RegionAllocatorException) if
    this $(D RegionAllocator) is not the most recently created still-existing
    $(D RegionAllocator) using its $(D RegionAllocatorStack) instance.
    */
    void freeLast() {
        mixin(getStackImplMixin);
        if(impl.regionIndex != this.correctRegionIndex) {
            throw new RegionAllocatorException(
                "Cannot free memory to a RegionAllocator that is not " ~
                "currently at the top of the stack, or memory that has not " ~
                "been allocated with this instance.",
                __FILE__,
                __LINE__
            );
        }

        with(*impl) {
            auto lastPos = bookkeep[--bookkeepIndex].p;

            // Handle large blocks.
            if(lastPos > space + segmentSize || lastPos < space) {
                alignedFree(lastPos);
                impl.nLargeObjects--;
                return;
            }

            used = (cast(size_t) lastPos) - (cast(size_t) space);
            if (nBlocks > 1 && used == 0) {
                impl.freeSegment();
            }
        }
    }

    /**
    Checks that $(D ptr) is a pointer to the block that would be freed by
    $(D freeLast) then calls $(D freeLast).  Throws a 
    $(D RegionAllocatorException) if the pointer does not point to the
    block that would be freed by $(D freeLast).
    */
    void free(void* ptr) {
        mixin(getStackImplMixin);
        auto lastPos = impl.bookkeep[impl.bookkeepIndex - 1].p;
        if(ptr !is lastPos) {
            throw new RegionAllocatorException(
                "Cannot free RegionAllocator memory in non-LIFO order.",
                __FILE__,
                __LINE__
            );
        }
        
        freeLast();
    }
    
    /**
    Attempts to resize a previously allocated block of memory in place.  
    This is possible only if $(D ptr) points to the beginning of the last
    block allocated by this $(D RegionAllocator) instance and, in the
    case where $(D newSize) is greater than the old size, there is
    additional space in the segment that $(D ptr) was allocated from.  
    Additionally, blocks larger than this $(D RegionAllocator)'s segment size
    cannot be grown or shrunk.
    
    Returns:  True if the block was successfully resized, false otherwise.
    */
    bool resize(const(void)* ptr, size_t newSize) {
        mixin(getStackImplMixin);
        
        // This works since we always allocate in increments of alignBytes
        // no matter what the allocation size.
        newSize = allocSize(newSize);
        
        with(*impl) {
            auto lastPos = bookkeep[bookkeepIndex - 1].p;
            if(ptr !is lastPos) {
                return false;
            }

            // If it's a large block directly allocated on the heap,
            // then we definitely can't resize it in place.
            if(lastPos > space + segmentSize || lastPos < space) {
                return false;
            }           
        
            immutable blockSize = used - ((cast(size_t) lastPos) -
                cast(size_t) space);
            immutable sizediff_t diff = newSize - blockSize;
            
            if(cast(sizediff_t) (impl.segmentSize - used) >= diff) {
                used += diff;
                return true;
            }
        }
        
        return false;
    }
    
    /**
    Returns whether the $(D RegionAllocatorStack) used by this
    $(D RegionAllocator) instance is scanned by the garbage collector.
    */
    bool gcScanned() @property const pure nothrow @safe {
        return stack.gcScanned;
    }

    /**Allocates an array of type $(D T).  $(D T) may be a multidimensional 
    array.  In this case sizes may be specified for any number of dimensions 
    from 1 to the number in $(D T).

    Examples:
    ---
    auto alloc = newRegionAllocator();
    double[] arr = alloc.newArray!(double[])(100);
    assert(arr.length == 100);

    double[][] matrix = alloc.newArray!(double[][])(42, 31);
    assert(matrix.length == 42);
    assert(matrix[0].length == 31);
    ---
    */
    auto newArray(T, I...)(I sizes)
    if(allSatisfy!(isIntegral, I)) {

        static void initialize(R)(R toInitialize) {
            static if(isArray!(ElementType!R)) {
                foreach(elem; toInitialize) {
                    initialize(elem);
                }
            } else {
                toInitialize[] = ElementType!(R).init;
            }
        }

        auto ret = uninitializedArray!(T, I)(sizes);
        initialize(ret);
        return ret;
    }

    /**
    Same as $(D newArray), except skips initialization of elements for
    performance reasons.
    */
    auto uninitializedArray(T, I...)(I sizes)
    if(allSatisfy!(isIntegral, I)) {
        static assert(sizes.length >= 1,
            "Cannot allocate an array without the size of at least the first " ~
            " dimension.");
        static assert(sizes.length <= nDimensions!T,
            to!string(sizes.length) ~ " dimensions specified for a " ~
            to!string(nDimensions!T) ~ " dimensional array.");

        alias typeof(T.init[0]) E;

        auto ptr = cast(E*) allocate(sizes[0] * E.sizeof);
        auto ret = ptr[0..sizes[0]];

        static if(sizes.length > 1) {
            foreach(ref elem; ret) {
                elem = uninitializedArray!(E)(sizes[1..$]);
            }
        }

        return ret;
    }

    /**
    Returns the number of bytes to which an allocation of size nBytes is
    guaranteed to be aligned.
    */
    static size_t alignBytes(size_t nBytes) {
        return .alignBytes;
    }

    /**
    Returns the number of bytes used to satisfy an allocation request
    of $(D nBytes).  Will return a value greater than or equal to
    $(D nBytes) to account for alignment overhead.
    */
    static size_t allocSize(size_t nBytes) pure nothrow {
        static assert(isPowerOf2(.alignBytes));
        return (nBytes + (.alignBytes - 1)) & (~(.alignBytes - 1));
    }

    /**
    False because memory allocated by this allocator is not automatically
    reclaimed by the garbage collector.
    */
    enum isAutomatic = false;

    /**
    True because, when the last last copy of a $(RegionAllocator) instance
    goes out of scope, the memory it references is automatically freed.
    */
    enum isScoped = true;

    /**
    True because if memory is freed via $(D free()) instead of $(D freeLast()) 
    then the pointer is checked for validity.
    */
    enum freeIsChecked = true;

    /**
    Returns the segment size of the $(D RegionAllocatorStack) used by this
    $(D RegionAllocator).
    */
    size_t segmentSize() @property {
        mixin(getStackImplMixin);
        return impl.segmentSize;
    }

    /**
    Returns the maximum number of bytes that may be allocated in the
    current segment.
    */
    size_t segmentSlack() @property {
        mixin(getStackImplMixin);
        return impl.segmentSize - impl.used;
    }

    /**
    Copies $(D range) to an array.  The array will be located on the
    $(D RegionAllocator) stack if any of the following conditions apply:

    1.  $(D std.traits.hasIndirections!(ElementType!R)) is false.

    2.  $(D R) is a builtin array.  In this case $(D range) maintains pointers
        to all elements at least until $(D array) returns, preventing the
        elements from being freed by the garbage collector.  A similar 
        assumption cannot be made for ranges other than builtin arrays.
        
    3.  The $(D RegionAllocatorStack) instance used by this 
        $(D RegionAllocator) is scanned by the garbage collector.

    If none of these conditions is met, the array is returned on the C heap
    and $(D GC.addRange) is called.  In either case, $(D RegionAllocator.free),
    $(D RegionAllocator.freeLast), or the last copy of this $(D RegionAllocator)
    instance going out of scope will free the array as if it had been
    allocated on the $(D RegionAllocator) stack.

    Rationale:  The most common reason to call $(D array) on a builtin array is 
                to modify its contents inside a function without affecting the
                caller's view.  In this case $(D range) is not modified and
                prevents the elements from being freed by the garbage
                collector.  Furthermore, if the copy returned does need
                to be scanned, the client can call $(D GC.addRange) before
                modifying the original array.

    Examples:
    ---
    auto alloc = newRegionAllocator();
    auto arr = alloc.array(iota(5));
    assert(arr == [0, 1, 2, 3, 4]);
    ---
    */
    Unqual!(ElementType!(R))[] array(R)(R range) if(isInputRange!R) {
        alias Unqual!(ElementType!(R)) E;
        if(gcScanned || !hasIndirections!E || isArray!R) {
            return arrayImplStack(range);
        } else {
            return arrayImplHeap(range);
        }
    }    
}

/**
Returns a new $(D RegionAllocator) that uses the default thread-local 
$(D RegionAllocatorStack) instance.
*/
RegionAllocator newRegionAllocator() {
    return RegionAllocator(getThreadLocal());
}

// Finishes copying a range to a C heap allocated array.  Assumes the first
// half of the input array is stuff already copied and the second half is
// free space.
private void finishCopy(T, U)(ref T[] result, U range, size_t alreadyCopied) {
    void doRealloc() {
        auto newPtr = cast(T*) alignedRealloc(
            result.ptr, result.length * T.sizeof * 2, result.length * T.sizeof
        );
        result = newPtr[0..result.length * 2];
    }

    auto index = alreadyCopied;
    foreach(elem; range) {
        if(index == result.length) doRealloc();
        result[index++] = elem;
    }

    result = result[0..index];
}

unittest {
    auto alloc = newRegionAllocator();
    auto arr = alloc.array(iota(5));
    assert(arr == [0, 1, 2, 3, 4]);

    // Create quick and dirty finite but lengthless range.
    static struct Count {
        uint num;
        uint upTo;
        @property size_t front() {
            return num;
        }
        void popFront() {
            num++;
        }
        @property bool empty() {
            return num >= upTo;
        }
    }

    alloc.allocate(1024 * 1024 * 3);
    Count count;
    count.upTo = 1024 * 1025;
    auto asArray = alloc.array(count);
    foreach(i, elem; asArray) {
        assert(i == elem, to!(string)(i) ~ "\t" ~ to!(string)(elem));
    }
    assert(asArray.length == 1024 * 1025);
    alloc.freeLast();
    alloc.freeLast();
    
    while(alloc.stack.impl.refCountedPayload.freeList.index > 0) {
        alignedFree(alloc.stack.impl.refCountedPayload.freeList.pop());
    }
}

unittest {
    auto alloc = newRegionAllocator();
    double[] arr = alloc.uninitializedArray!(double[])(100);
    assert(arr.length == 100);

    double[][] matrix = alloc.uninitializedArray!(double[][])(42, 31);
    assert(matrix.length == 42);
    assert(matrix[0].length == 31);

    double[][] mat2 = alloc.newArray!(double[][])(3, 1);
    assert(mat2.length == 3);
    assert(mat2[0].length == 1);

    import std.math;
    assert(isNaN(mat2[0][0]));
}

unittest {
    /* Not a particularly good unittest in that it depends on knowing the
     * internals of RegionAllocator, but it's the best I could come up w/.  This
     * is really more of a stress test/sanity check than a normal unittest.*/

    // Make sure state is completely reset.
    clear(threadLocalStack.impl);
    threadLocalStack = RegionAllocatorStack.init;
    threadLocalInitialized = false;

     // First test to make sure a large number of allocations does what it's
     // supposed to in terms of reallocing bookkeep[], etc.
     enum nIter =  defaultSegmentSize * 5 / alignBytes;

    {
         auto alloc = newRegionAllocator();
         foreach(i; 0..nIter) {
             alloc.allocate(alignBytes);
         }
         assert(alloc.stack.impl.refCountedPayload.nBlocks == 5,
            to!string(alloc.stack.impl.refCountedPayload.nBlocks));
         assert(alloc.stack.impl.refCountedPayload.nFree == 0);
         foreach(i; 0..nIter) {
            alloc.freeLast();
        }
        assert(alloc.stack.impl.refCountedPayload.nBlocks == 1);
        assert(alloc.stack.impl.refCountedPayload.nFree == 2);

        // Make sure logic for freeing excess blocks works.  If it doesn't this
        // test will run out of memory.
        enum allocSize = defaultSegmentSize / 2;
        foreach(i; 0..50) {
            foreach(j; 0..50) {
                alloc.allocate(allocSize);
            }
            foreach(j; 0..50) {
                alloc.freeLast();
            }
        }

        // Make sure data is stored properly.
        foreach(i; 0..10) {
            alloc.allocate(allocSize);
        }
        foreach(i; 0..5) {
            alloc.freeLast();
        }
        void* space = alloc.stack.impl.refCountedPayload.space;
        size_t used = alloc.stack.impl.refCountedPayload.used;

        {
            auto alloc2 = newRegionAllocator();
            auto arrays = new uint[][10];
            
            foreach(i; 0..10) {
                uint[] data = alloc2.uninitializedArray!(uint[])(250_000);
                foreach(j, ref e; data) {
                    e = cast(uint) (j * (i + 1));  
                                    }
                arrays[i] = data;
            }

            // Make stuff get overwrriten if blocks are getting GC'd when
            // they're not supposed to.
            GC.minimize;  // Free up all excess pools.
            uint[][] foo;
            foreach(i; 0..40) {
                foo ~= new uint[1_048_576];
            }
            foo = null;

            for(size_t i = 9; i != size_t.max; i--) {
                foreach(j, e; arrays[i]) {
                    assert(e == j * (i + 1));
                }
            }
        }

        assert(space == alloc.stack.impl.refCountedPayload.space,
            text(space, '\t', alloc.stack.impl.refCountedPayload.space));
        assert(used == alloc.stack.impl.refCountedPayload.used,
            text(used, '\t', alloc.stack.impl.refCountedPayload.used));
        while(alloc.stack.impl.refCountedPayload.nBlocks > 1 || 
        alloc.stack.impl.refCountedPayload.used > 0) {
            alloc.freeLast();
        }
    }

    // Test that everything is really getting destroyed properly when
    // destroy() is called.  If not then this test will run out of memory.
    foreach(i; 0..1000) {
        clear(threadLocalStack.impl);
        threadLocalInitialized = false;

        auto alloc = newRegionAllocator();
        foreach(j; 0..1_000) {
            auto ptr = alloc.allocate(20_000);
            assert((cast(size_t) ptr) % alignBytes == 0);
        }

        foreach(j; 0..500) {
            alloc.freeLast();
        }
    }
}

unittest {
    // Make sure the basics of using explicit stacks work.
    auto stack = RegionAllocatorStack(4 * 1024 * 1024, GCScan.no);
    auto alloc = stack.newRegionAllocator();
    auto arr = alloc.array(iota(5));
    assert(arr == [0, 1, 2, 3, 4]);
    auto ptr = alloc.allocate(5);
    
    auto alloc2 = newRegionAllocator();
    auto ptr2 = alloc2.allocate(5);
    auto ptr3 = alloc.allocate(5);
}

unittest {
    // Make sure the stacks get freed properly when they go out of scope.
    // If they don't then this will run out of memory.
    foreach(i; 0..100_000) {
        auto stack = RegionAllocatorStack(4 * 1024 * 1024, GCScan.no);
    }
}

unittest {
    // Make sure resize works properly.
    auto alloc = newRegionAllocator();
    auto arr1 = alloc.array(iota(4));
    auto res = alloc.resize(arr1.ptr, 8 * int.sizeof);
    assert(res);
    arr1 = arr1.ptr[0..8];
    copy(iota(4, 8), arr1[4..$]);
    
    auto arr2 = alloc.newArray!(int[])(8);
    
    // If resizing resizes to something too small, this will have been 
    // overwritten:
    assert(arr1 == [0, 1, 2, 3, 4, 5, 6, 7], text(arr1));
    
    alloc.free(arr2.ptr);
    auto res2 = alloc.resize(arr1.ptr, 4 * int.sizeof);
    assert(res2);
    arr2 = alloc.newArray!(int[])(8);
    
    // This time the memory in arr1 should have been overwritten.
    assert(arr1 == [0, 1, 2, 3, 0, 0, 0, 0]);
}

unittest {
    // Make sure that default thread-local stacks get freed properly at the
    // termination of a thread.  If they don't then this will run out of
    // memory.
    
    import core.thread;
    foreach(i; 0..100) {
        auto t = new Thread({ 
            threadLocalSegmentSize = 100 * 1024 * 1024;
            newRegionAllocator(); 
        });
        t.start();
        t.join();
    }
}

unittest {
    // Make sure assignment works as advertised.
    RegionAllocator alloc;
    auto alloc2 = newRegionAllocator();
    auto ptr = alloc2.allocate(8);
    alloc = alloc2;
    alloc.freeLast();
    auto ptr2= alloc2.allocate(8);
    assert(ptr is ptr2);
}

 // Simple, fast stack w/o error checking.
static struct SimpleStack(T) { 
    private size_t capacity;
    private size_t index;
    private T* data;
    private enum sz = T.sizeof;

    private static size_t max(size_t lhs, size_t rhs) pure nothrow {
        return (rhs > lhs) ? rhs : lhs;
    }

    void push(T elem) {
        if (capacity == index) {
            capacity = max(16, capacity * 2);
            data = cast(T*) core.stdc.stdlib.realloc(data, capacity * sz);
        }
        data[index++] = elem;
    }

    T pop() {
        index--;
        auto ret = data[index];
        return ret;
    }

    void destroy() {
        if(data) {
            core.stdc.stdlib.free(data);
            data = null;
        }
    }
}

private  void outOfMemory()  {
    throw new OutOfMemoryError("Out of memory in RegionAllocator.");
}

// Memory allocation routines.  These wrap allocate(), free() and realloc(),
// and guarantee alignment.
private enum size_t alignBytes = 16;

private void* alignedMalloc(size_t size, bool shouldAddRange = false) {
    // We need (alignBytes - 1) extra bytes to guarantee alignment, 1 byte
    // to store the shouldAddRange flag, and ptrSize bytes to store
    // the pointer to the beginning of the block.
    void* toFree = core.stdc.stdlib.malloc(
        alignBytes + ptrSize + size
    );

    if(toFree is null) outOfMemory();

    // Add the offset for the flag and the base pointer.
    auto intPtr = cast(size_t) toFree + ptrSize + 1;

    // Align it.
    intPtr = (intPtr + alignBytes - 1) & (~(alignBytes - 1));
    auto ret = cast(void**) intPtr;

    // Store base pointer.
    (cast(void**) ret)[-1] = toFree;

    // Store flag.
    (cast(bool*) ret)[-1 - ptrSize] = shouldAddRange;

    if(shouldAddRange) {
        GC.addRange(ret, size);
    }

    return ret;
}

private void alignedFree(void* ptr) {
    // If it was allocated with alignedMalloc() then the pointer to the
    // beginning is at ptr[-1].
    auto addedRange = (cast(bool*) ptr)[-1 - ptrSize];

    if(addedRange) {
        GC.removeRange(ptr);
    }

    core.stdc.stdlib.free( (cast(void**) ptr)[-1]);
}

// This is used by RegionAllocator, but I'm not sure enough that its interface
// isn't going to change to make it public and document it.
private void* alignedRealloc(void* ptr, size_t newLen, size_t oldLen) {
    auto storedRange = (cast(bool*) ptr)[-1 - ptrSize];
    auto newPtr = alignedMalloc(newLen, storedRange);
    memcpy(newPtr, ptr, oldLen);

    alignedFree(ptr);
    return newPtr;
}

private union SizetPtr {
    size_t i;
    void* p;
}
