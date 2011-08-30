/** Provides SplineView structure.
  *
  * Version: 0.4-a
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev.
  * License: $(WEB boost.org/LICENSE_1_0.txt, Boost License 1.0)
  */
module scid.splines.univariate.view;

import scid.common.traits;
import scid.splines.support;

/** Structure to access a one-dimensional spline, which provides the interface
  * for function evaluation and accelerates this action if multiple successive
  * evaluations are made on the same interval.
  *
  * Multiple SplineView structures can refer the same spline.
  *
  * Under the hood, it caches the previous value of an index lookup and return
  * it immediately if the subsequent interpolation point falls in the same
  * interval.
  *
  * Params:
  *     Tspline = type of spline to wrap
  */
struct SplineView(Tspline)
{
    // Pointer to the wrapped spline
    private Tspline* _spline; // TODO: make it a property?

    // Types from the wrapped spline
    alias BaseElementType!(typeof(_spline._f)) FuncType;
    alias BaseElementType!(typeof(_spline._x)) VarType;

    // Index lookup
    private
    {
        size_t _index = 0; // The index of current interval

        // Find the value x in a sorted array a
        void _updateIndex(VarType x)
        {
            // Does current interval contain x?
            if((_spline._x[_index + 1] < x) || (_spline._x[_index] > x))
                // Find new interval if it doesn't
                _index = binarySearch(_spline._x, x);
        }
    }

    public
    {
        /** Wrap the specified spline
          */
        this(ref Tspline spline)
        {
            _spline = &spline;
        }

        /** Evaluate function value at a given point
          * Params:
          *     x = variable value
          */
        FuncType eval(VarType x)
        in
        {
            assert(_spline.pointInsideRange(x),
                   "Variable value out of interpolation range");
        }
        body
        {
            _updateIndex(x);
            return _spline._calcFunction(x, _index);
        }

        /** Evaluate the first derivative at a given point
          * Params:
          *     x = variable value
          */
        FuncType deriv(VarType x)
        in
        {
            assert(_spline.pointInsideRange(x),
                   "Variable value out of interpolation range");
        }
        body
        {
            _updateIndex(x);
            return _spline._calcDeriv(x, _index);
        }

        // TODO: derivatives and integration
    }
}
