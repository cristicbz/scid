/** Contains features that are used by all splines but need not to be visible
  * from outside the library.
  *
  * Warning: This module is for scid.splines internal use only.
  * Its interface can be changed at any time and no external program should
  * rely on it.
  *
  * Version: 0.7-b
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev.
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  */
module scid.splines.support;

public import scid.common.traits;
public import scid.internal.regionallocator;

/** Returns array type based on given type.
  * If given type is already array-like does nothing
  */
template ProduceArray(ElementOrContainer)
{
    static if(is(BaseElementType!ElementOrContainer == ElementOrContainer))
        alias ElementOrContainer[] ProduceArray;
    else
        alias ElementOrContainer ProduceArray;
}

unittest
{
    static assert(is(ProduceArray!double == double[]));
    static assert(is(ProduceArray!(double[]) == double[]));
}

/** Common binary search in a sorted array.
  * Returns index of the first element in the first interval that contains x
  */
size_t binarySearch(A, X)(A a, X x)
{
    size_t ilo = 0;
    size_t ihi = a.length;
    while(ihi - ilo > 1)
    {
        size_t i = (ihi + ilo) / 2;
        if(a[i] == x)
        {
            // Exact match
            ilo = i;
            ihi = ilo + 1; // This is to violate loop condition
        }
        else if(a[i] > x)
            ihi = i; // x is somwhere on the right
        else
            ilo = i; // x is somwhere on the left
    }
    if(ilo == a.length - 1)
        // if x is equal to the last element return the last interval
        --ilo;
    return ilo;
}

unittest
{
    double[] a = [0, 1, 2, 3];
    assert(binarySearch(a, 0) == 0);
    assert(binarySearch(a, 0.5) == 0);
    assert(binarySearch(a, 1.5) == 1);
    assert(binarySearch(a, 2.5) == 2);
    assert(binarySearch(a, 3) == 2);
}

/* Create a RegionAllocator using the stack specified.
 * If null is specified then thread-local stack is used.
 */
RegionAllocator newRegionAllocatorInStack(RegionAllocatorStack* stack)
{
    if(stack)
        return stack.newRegionAllocator();
    else
        return newRegionAllocator();
}
