/** Contains common features of boundary conditions of univariate splines.
  *
  * Version: 0.6-a
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev.
  * License: $(WEB boost.org/LICENSE_1_0.txt, Boost License 1.0)
  */
module scid.splines.univariate.boundcond;

/// Boundary conditons (BC) types
enum BoundCond
{
    /** Natural boundary conditions. BC values are ignored.
      * See particular splines' descriptions for details.
      */
    natural,
    /** Periodic spline. BC values are ignored.
      * Setting it on one side automatically sets it on the other.
      */
    periodic,
    /** Corresponding BC value is the first derivative */
    deriv1,
    /** Corresponding BC value is the second derivative */
    deriv2
}

package mixin template boundaryConditions()
{
    private BoundCond _bcLeftType;

    /// Type of BC on the left side
    @property void bcLeftType(BoundCond bc)
    in
    {
        assert(bcIsSupported(bc));
    }
    body
    {
        _bcLeftType = bc;
        static if(bcIsSupported(BoundCond.periodic))
        {
            if(bc == BoundCond.periodic)
                _bcRightType = BoundCond.periodic;
        }
    }

    /// ditto
    @property BoundCond bcLeftType()
    {
        return _bcLeftType;
    }

    private FuncType _bcLeftVal;

    /// Value of BC on the left side
    @property void bcLeftVal(FuncType val)
    {
        _bcLeftVal = val;
    }

    /// ditto
    @property FuncType bcLeftVal()
    {
        return _bcLeftVal;
    }

    private BoundCond _bcRightType;

    /// Type of BC on the right side
    @property void bcRightType(BoundCond bc)
    in
    {
        assert(bcIsSupported(bc));
    }
    body
    {
        _bcRightType = bc;
        static if(bcIsSupported(BoundCond.periodic))
        {
            if(bc == BoundCond.periodic)
                _bcLeftType = BoundCond.periodic;
        }
    }

    /// ditto
    @property BoundCond bcRightType()
    {
        return _bcRightType;
    }

    private FuncType _bcRightVal;

    /// Value of BC on the right side
    @property void bcRightVal(FuncType val)
    {
        _bcRightVal = val;
    }

    /// ditto
    @property FuncType bcRightVal()
    {
        return _bcRightVal;
    }
}
