/** Provides univariate Akima spline.
  *
  * Version: 0.7-b
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev.
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  * Bugs:
  *     No optimization support. $(BR)
  *     1st and 2nd derivative boundary conditions are not supported.
  */
module scid.splines.univariate.akima;

public import scid.splines.univariate.boundcond;

import std.math;

import scid.common.meta;
import scid.splines.univariate.base;
import scid.splines.univariate.poly3;

/** One-dimensional Akima interpolation.
  *
  * Natural boundary conditions are non-trivial and described in the original
  * paper [H.Akima, Journal of the ACM, 17(4), 589 (1970)]
  *
  * Params:
  *     EocVar = type of variable or variable value array (VVA)
  *     EocFunc = type of function or function value array (FVA)
  *     optim = spline optimization type (normal by default)
  */
struct SplineAkima(EocVar, EocFunc,
                   SplineOptim optim = SplineOptim.normal)
{
    // FIXME: Add different boundary conditions
    // TODO: Add a mechanism for adding points to the curve
    // FIXME: Implement support of compound types
    // TODO: Unittest

    mixin splineBase!(EocVar, EocFunc);
    mixin poly3Base!(VarType, FuncType);

    // Interpolant storage
    private
    {
        static if(optim == SplineOptim.fixVar)
        {
            /* Data depending only on variable values */
            // FIXME: implement after testing of the algorithm

            void _allocOptFixVar(size_t maxSize)
            {
            }
        }
    }

    // Numereical core
    package
    {
        static if(optim == SplineOptim.fixVar)
        {
            // TODO: implement this area when the algorithm will be tested
            static assert(false, "fixVar mode is not implemented yet");

            void _calcVarDependent()
            {
                // FIXME: implement after testing of the algorithm
            }

            void _calcFuncDependent()
            {
                // FIXME: implement after testing of the algorithm
            }
        }
        else
        {
            void _calcAll()
            {
                calcCoeffs(_x, _f,
                           _bcLeftType, _bcRightType,
                           _bcLeftVal, _bcRightVal,
                           _c1, _c2, _c3,
                           workspaceRegAllocStack);
            }
        }
    }

    // Boundary conditions support
    public
    {
        static bool bcIsSupported(BoundCond bc)
        {
            switch(bc)
            {
                case BoundCond.natural: return true;
                case BoundCond.periodic: return true;
                case BoundCond.deriv1: return true;
                case BoundCond.deriv2: return false;
                default: return false;
            }
        }

        mixin boundaryConditions;
    }

    public
    {
        /// Minimal number of points needed for the spline
        enum size_t minPoints = 3;

        /** Reserve memory for the spline.
          */
        this(size_t maxSize)
        {
            _allocContents(maxSize);
        }

        /** Create spline for given variable and function value arrays.
          *
          * Params:
          *     x = variable value array
          *     y = function value array
          *     calcNow = whether to calculate the spline immediately
          *               (true by default)
          */
        this(VarArray x, FuncArray y, bool calcNow = true,
             BoundCond bcLType = BoundCond.natural,
             BoundCond bcRType = BoundCond.natural,
             FuncType bcLVal = Zero!(FuncType),
             FuncType bcRVal = Zero!(FuncType))
        {
            this(x.length);
            setAll(x, y);
            /* For the following case: this(x, y, calcNow, periodic) */
            if(bcLType == BoundCond.periodic)
                bcRType = BoundCond.periodic;
            _bcLeftType = bcLType;
            _bcRightType = bcRType;
            _bcLeftVal = bcLVal;
            _bcRightVal = bcRVal;
            if(calcNow)
                calculate();
        }
    }
}

/* Calculate Akima spline coefficients
 * [H.Akima, Journal of the ACM, 17(4), 589 (1970)]
 */
private void calcCoeffs(VarType, FuncType)
                       (in VarType[] x,
                        in FuncType[] f,
                        in BoundCond bcLeftType,
                        in BoundCond bcRightType,
                        in FuncType bcLeftVal,
                        in FuncType bcRightVal,
                        FuncType[] c1,
                        FuncType[] c2,
                        FuncType[] c3,
                        RegionAllocatorStack* wsras)
{
    // The index of the last point of the spline
    size_t N = x.length - 1;
    auto workspace = newRegionAllocatorInStack(wsras);
    auto d = workspace.uninitializedArray!(FuncType[])(N + 2);
    auto w = workspace.uninitializedArray!(FuncType[])(N + 3);

    // Calculate slopes
    for(size_t i = 1; i <= N; ++i)
        d[i] = (f[i] - f[i - 1]) / (x[i] - x[i - 1]);
    /* Process boundary points:
     *   calculate in a special way d[0] and d[N + 1] slopes
     *   for which x[-1] and x[N + 1] would have been required
     */
    if(bcLeftType == BoundCond.periodic)
    {
        d[0] = d[N];
        d[N + 1] = d[1];
    }
    else
    {
        if(bcLeftType == BoundCond.natural)
            d[0] = 2 * d[1] - d[2];
        else if(bcLeftType == BoundCond.deriv1)
            d[0] = 2 * bcLeftVal - d[1];

        if(bcRightType == BoundCond.natural)
            d[N + 1] = 2 * d[N] - d[N - 1];
        else if(bcRightType == BoundCond.deriv1)
            d[N + 1] = 2 * bcRightVal - d[N];
    }

    // Calculate weights
    for(size_t i = 1; i <= N + 1; ++i)
        w[i] = abs(d[i] - d[i - 1]);
    /* Process boundary points:
     *   calculate in a special way w[0] and w[N + 2] weights
     *   for which d[-1] and d[N + 2] would have been required
     */
    if(bcLeftType == BoundCond.periodic)
    {
        w[0] = w[N];
        w[N + 2] = w[2];
    }
    else
    {
        if(bcLeftType == BoundCond.natural)
            w[0] = w[1];
        else if(bcLeftType == BoundCond.deriv1)
            w[0] = w[2];

        if(bcRightType == BoundCond.natural)
            w[N + 2] = w[N + 1];
        else if(bcRightType == BoundCond.deriv1)
            w[N + 2] = w[N];
    }
    // Calculate the first derivatives
    for(size_t i = 0; i <= N; ++i)
        if(w[i] + w[i + 2] > 0)
            c1[i] = (w[i] * d[i] + w[i + 2] * d[i + 1])
                     / (w[i] + w[i + 2]);
        else
            c1[i] = (d[i] + d[i + 1]) / 2;
    // Calculate the remaining coefficients
    for(size_t i = 0; i < N; ++i)
    {
        VarType h = x[i + 1] - x[i];
        VarType v = f[i + 1] - f[i];
        c2[i] = (3 * v - h * (2 * c1[i] + c1[i + 1]))
                 / (h * h);
        c3[i] = (-2 * v + h * (c1[i] + c1[i + 1]))
                 / (h * h * h);
    }
}
