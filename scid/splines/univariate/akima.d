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

/** One-dimensional Akima interpolation.
  *
  * Natural boundary conditions are non-trivial and described in
  * [H.Akima, Journal of the ACM, 17(4), 589 (1970)]
  *
  * Params:
  *     EocVar = type of variable or variable value array (VVA)
  *     EocFunc = type of function or function value array (FVA)
  *     optim = spline optimization type (normal by default)
  */
struct SplineAkima(EocVar, EocFunc,
                   SplineOptim optim = SplineOptim.normal)
{
    // TODO: Add different boundary conditions
    // TODO: Add a mechanism for adding points to the curve
    // TODO: Implement support of compound types
    // TODO: Unittest

    mixin splineBase!(EocVar, EocFunc);

    // Data and workspaces
    private
    {
        // Spline parameters:
        FuncType[] _c1;
        FuncType[] _c2;
        FuncType[] _c3;
        /* The interpolant is:
         *     f(x) = _f[i] + _c1[i] * dx + _c2[i] * dx*dx + _c3[i] * dx*dx*dx
         *     dx = x - _x[i]
         */

        void _allocContents(size_t maxSize) // TODO: use scid.core.memory
        {
            _c1.length = maxSize;
            _c2.length = maxSize - 1;
            _c3.length = maxSize - 1;

            static if(optim == SplineOptim.fixVar)
            {
                // FIXME: implement after testing of the algorithm
            }
            _maxSize = maxSize;
        }
    }

    // Numereical core
    package
    {
        static if(optim == SplineOptim.fixVar)
        {
            // TODO: implement this area when the algorithm will be tested
            static assert(false, "fixVar mode is not implemented yet");

            bool _needUpdateVar;

            /* Data depending only on variable values */
            // FIXME: implement after testing of the algorithm

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

    // Function evaluation code
    package
    {
        FuncType _calcFunction(VarType x, size_t index)
        {
            double dx = x - _x[index];
            return _f[index] + dx * (_c1[index]
                                     + dx * (_c2[index] + dx * _c3[index]));
        }

        // Calculate first derivative in a given interval
        FuncType _calcDeriv(VarType x, size_t index)
        {
            double dx = x - _x[index];
            return _c1[index] + dx * (2 * _c2[index] + dx * 3 * _c3[index]);
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
                case BoundCond.deriv1: return false;
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

        /** Spline coefficients (read only)
          */
        @property const(FuncType)[] c1()
        {
            return _c1[0..(pointsNum - 1)];
        }

        ///ditto
        @property const(FuncType)[] c2()
        {
            return _c2[0..(pointsNum - 1)];
        }

        ///ditto
        @property const(FuncType)[] c3()
        {
            return _c3[0..(pointsNum - 1)];
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
    // FIXME: real only
    // The index of the last point of the spline
    size_t N = x.length - 1;
    auto workspace = newRegionAllocatorInStack(wsras);
    auto d = workspace.uninitializedArray!(FuncType[])(N + 2);
    auto w = workspace.uninitializedArray!(FuncType[])(N + 3);

    // Calculate slopes
    for(size_t i = 1; i <= N; ++i)
        d[i] = (f[i] - f[i - 1]) / (x[i] - x[i - 1]);
    // Process boundary points
    if(bcLeftType == BoundCond.periodic)
    {
        d[0] = d[N];
        d[N + 1] = d[1];
    }
    else
    {
        if(bcLeftType == BoundCond.natural)
            d[0] = 2 * d[1] - d[2];

        if(bcRightType == BoundCond.natural)
            d[N + 1] = 2 * d[N] - d[N - 1];
    }

    // Calculate weights
    for(size_t i = 1; i <= N + 1; ++i)
        w[i] = abs(d[i] - d[i - 1]);
    // Process boundary points
    if(bcLeftType == BoundCond.periodic)
    {
        w[0] = w[N];
        w[N + 2] = w[2];
    }
    else
    {
        if(bcLeftType == BoundCond.natural)
            w[0] = w[1];

        if(bcRightType == BoundCond.natural)
            w[N + 2] = w[N + 1];
    }
    w[0] = w[1];
    w[N + 2] = w[N + 1];
    // Calculate the first derivatives
    for(size_t i = 0; i <= N; ++i)
        if(w[i] + w[i + 2] > 0)
            c1[i] = (w[i] * d[i] + w[i + 2] * d[i + 1])
                     / (w[i] + w[i + 2]);
        else
            c1[i] = (d[i] + d[i + 2]) / 2;
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
