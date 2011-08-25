module scid.splines.base;

public import scid.common.traits;

/// The type of optimization
enum SplineOptim
{
    normal, /// No special features.
    fixVar /** Preliminary calculations depending only on VVA are made if it
             * changes. Minimal amount of operations is performed when FVA
             * changes.
             */
}

/* Returns array type based on given type.
   If given type is already array-like does nothing */
/* NOTE: visibility levels do not work for templates!
 */
package template ProduceArray(ElementOrContainer)
{
    static if(is(BaseElementType!ElementOrContainer == ElementOrContainer))
        alias ElementOrContainer[] ProduceArray;
    else
        alias ElementOrContainer ProduceArray;
}

/* Common binary search in a sorted array.
   Returns index of the first element in the first interval that contains x */
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
