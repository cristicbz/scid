/** Provides univariate cubic spline.
  *
  * Version: 0.5-a
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev.
  * License: $(WEB boost.org/LICENSE_1_0.txt, Boost License 1.0)
  */
module scid.splines.univariate.cubic;

import scid.common.meta;
import scid.splines.univariate.base;

/// Boundary conditons (BC) types
enum BoundCond
{
    /** Periodic spline. BC values are ignored.
      * Setting it on one side automatically sets it on the other.
      */
    periodic,
    /** Corresponding BC value is the first derivative */
    deriv1,
    /** Corresponding BC value is the second derivative */
    deriv2
}

/** Cubic one-dimensional spline (order = 3, defect = 1).
  *
  * Params:
  *     EocVar = type of variable
  *     EocFunc = type of function
  *     optim = spline optimization type (normal by default)
  */
struct SplineCubic(EocVar, EocFunc,
                   SplineOptim optim = SplineOptim.normal)
{
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
            _c1.length = maxSize - 1;
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
            // Calculate the interpolant parameters
            void _calcAll()
            {
                if(_bcLeftType == BoundCond.periodic)
                {
                    calcCoeffsPeriodic(_x, _f,
                                       _c1, _c2, _c3,
                                       workspaceRegAllocStack);
                }
                else
                {
                    calcCoeffs(_x, _f,
                               _bcLeftType, _bcRightType,
                               _bcLeftVal, _bcRightVal,
                               _c1, _c2, _c3,
                               workspaceRegAllocStack);
                }
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

    // Boundary conditions
    public
    {
        private BoundCond _bcLeftType = BoundCond.deriv2;

        void bcLeftType(BoundCond bc)
        {
            _bcLeftType = bc;
            if(bc == BoundCond.periodic)
                _bcRightType = BoundCond.periodic;
        }

        /// Type of BC on the left side
        BoundCond bcLeftType()
        {
            return _bcLeftType;
        }

        private FuncType _bcLeftVal = Zero!(FuncType);

        void bcLeftVal(FuncType val)
        {
            _bcLeftVal = val;
        }

        /// Value of BC on the left side
        FuncType bcLeftVal()
        {
            return _bcLeftVal;
        }

        private BoundCond _bcRightType = BoundCond.deriv2;

        void bcRightType(BoundCond bc)
        {
            _bcRightType = bc;
            if(bc == BoundCond.periodic)
                _bcLeftType = BoundCond.periodic;
        }

        /// Type of BC on the right side
        BoundCond bcRightType()
        {
            return _bcRightType;
        }

        private FuncType _bcRightVal = Zero!(FuncType);

        void bcRightVal(FuncType val)
        {
            _bcRightVal = val;
        }

        /// Value of BC on the right side
        FuncType bcRightVal()
        {
            return _bcRightVal;
        }
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
             BoundCond bcLType = BoundCond.deriv2,
             BoundCond bcRType = BoundCond.deriv2,
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

unittest
{
    alias SplineCubic!(double, double, SplineOptim.normal) MySpline;
    double[] x = [0, 1, 3];
    double[] y = [0, 1, 27];
    auto spl = MySpline(x, y, true,
                        BoundCond.deriv1,
                        BoundCond.deriv1,
                        0, 27);
    assert(spl._calcFunction(2, 1) == 8);
    assert(spl._calcDeriv(2, 1) == 12);
}

// Calculate cubic spline coefficients
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
in
{
    assert(f.length == x.length);
    assert(bcLeftType != BoundCond.periodic);
    assert(bcRightType != BoundCond.periodic);
    assert(c1.length >= x.length - 1);
    assert(c2.length >= x.length - 1);
    assert(c3.length >= x.length - 1);
}
body
{
    /* This function builds linear equation system
     * for coefficients c1[i]:
     *   a0[i] * c1[i-1] + a1[i] * c1[i] + a2[i] * c1[i+1] = b[i]
     * and then performs standard sweep procedure for
     * tridiagonal matrix.
     *
     * The coefficients and the right part of the system are
     * calculated on-the-fly
     */

    // The index of the last point of the spline
    size_t N = x.length - 1;
    // Allocate workspaces
    auto workspace = newRegionAllocatorInStack(wsras);
    auto wa = workspace.uninitializedArray!(VarType[])(N + 1);
    auto wb = workspace.uninitializedArray!(FuncType[])(N + 1);

    // Some variables for optimization
    VarType h0 = x[1] - x[0];
    FuncType v0 = f[1] - f[0];
    FuncType d0 = v0 / h0;
    // Linear equation coefficients
    VarType a0, a1, a2;
    // Linear equation right part
    FuncType b;

    // Boundary conditions on the left side (the first equation)
    if(bcLeftType == BoundCond.deriv1)
    {// FIXME: check
        a1 = 1;
        a2 = 0;
        b = bcLeftVal;
    }
    else if(bcLeftType == BoundCond.deriv2)
    {// FIXME: check
        a1 = 2;
        a2 = 1;
        b = 3 * d0 - h0 * bcLeftVal / 2;
    }
    // a0 won't be used in the first equation

    // Build'n'sweep forward

    /* The first equation (i=0) is processed separately since
     * there is no c1[-1] coefficient.
     */
    VarType factor = 1 / a1;
    wa[0] = -a2 * factor;
    wb[0] = b * factor;

    // Process non-border equations
    for(size_t i = 1; i < N; ++i)
    {
        // Optimization
        VarType h1 = x[i+1] - x[i];
        FuncType v1 = f[i+1] - f[i];
        FuncType d1 = v1 / h1;
        // Calculate eqaution coefficients
        a0 = h1 / (h0 + h1);
        a1 = 2;
        a2 = 1 - a0;
        b = 3 * (a0 * d0 + a2 * d1);
        // Forward sweep step
        factor = 1 / (a1 + a0 * wa[i - 1]);
        wa[i] = -a2 * factor;
        wb[i] = (b - a0 * wb[i - 1]) * factor;
        // Optimization
        h0 = h1;
        v0 = v1;
        d0 = d1;
    }

    // Boundary conditions on the right side (the last equation)
    if(bcRightType == BoundCond.deriv1)
    {// FIXME: check
        a0 = 0;
        a1 = 1;
        b = bcRightVal;
    }
    else if(bcRightType == BoundCond.deriv2)
    {// FIXME: check
        a0 = 1;
        a1 = 2;
        b = 3 * d0 - h0 * bcRightVal / 2;
    }

    // Sweep backward and calculate the coefficients

    /* The last equation (i=N) is dealt separately since
     * there is no c1[N] and c1[N+1] coefficients.
     */
    FuncType c1tmp = (b - a0 * wb[N - 1]) /
                     (a1 + a0 * wa[N - 1]);

    // Process non-border equations
    for(size_t i = N; i > 0; --i)
    {
        // Backward sweep step
        c1[i - 1] = wa[i - 1] * c1tmp + wb[i - 1];
        // Optimization
        VarType h = x[i] - x[i - 1];
        FuncType v = f[i] - f[i - 1];
        // Calculate remaining coefficients
        c2[i - 1] = (3 * v - h * (2 * c1[i - 1] + c1tmp))
                     / (h * h);
        c3[i - 1] = (-2 * v + h * (c1[i - 1] + c1tmp))
                     / (h * h * h);
        // Optimization
        c1tmp = c1[i - 1];
    }
}

/* Calculate periodic cubic spline coefficients
 * The last point in f[] is ignored and assumed to be equal to f[0]
 */
private void calcCoeffsPeriodic(VarType, FuncType)
                               (in VarType[] x,
                                in FuncType[] f,
                                FuncType[] c1,
                                FuncType[] c2,
                                FuncType[] c3,
                                RegionAllocatorStack* wsras)
in
{
    assert(f.length == x.length);
    assert(c1.length >= x.length - 1);
    assert(c2.length >= x.length - 1);
    assert(c3.length >= x.length - 1);
}
body
{
    /* This function builds linear equation system
     * for coefficients c1[i]:
     *   a0[i] * c1[i-1] + a1[i] * c1[i] + a2[i] * c1[i+1] = b[i]
     * and then performs standard cyclic sweep procedure for
     * tridiagonal matrix.
     *
     * The coefficients and the right part of the system are
     * calculated on-the-fly
     */

    // The index of the last point of the spline
    size_t N = x.length - 1;
    // Allocate workspaces
    auto workspace = newRegionAllocatorInStack(wsras);
    auto wa = workspace.uninitializedArray!(VarType[])(N + 1);
    auto wb = workspace.uninitializedArray!(FuncType[])(N + 1);
    auto wc = workspace.uninitializedArray!(VarType[])(N + 1);
    auto wu = workspace.uninitializedArray!(FuncType[])(N + 1);
    auto wv = workspace.uninitializedArray!(VarType[])(N + 1);

    // Some variables for optimization
    VarType h0 = x[1] - x[0];
    FuncType v0 = f[1] - f[0];
    FuncType d0 = v0 / h0;
    // Linear equation coefficients
    VarType a0, a1, a2;
    // Linear equation right part
    FuncType b;

    // Build'n'sweep forward

    /* The first equation (i=0) is dealt separately since
     * there is no _c1[-1] coefficient.
     */
    wa[0] = 0;
    wb[0] = Zero!(FuncType);
    wc[0] = 1;

    // OK, let's go!
    for(size_t i = 1; i < N; ++i)
    {
        // Optimization
        VarType h1 = x[i+1] - x[i];
        FuncType v1 = f[i+1] - f[i];
        FuncType d1 = v1 / h1;
        // Calculate eqaution coefficients
        a0 = h1 / (h0 + h1);
        a1 = 2;
        a2 = 1 - a0;
        b = 3 * (a0 * d0 + a2 * d1);
        // Forward sweep step
        VarType factor = 1 / (a1 + a0 * wa[i - 1]);
        wa[i] = -a2 * factor;
        wb[i] = (b - a0 * wb[i - 1]) * factor;
        wc[i] = -a0 * wc[i - 1] * factor;
        // Optimization
        h0 = h1;
        v0 = v1;
        d0 = d1;
    }

    // Sweep backward
    wu[N - 1] = wb[N - 1];
    wv[N - 1] = wa[N - 1] + wc[N - 1];
    for(size_t i = N - 1; i > 0; --i)
    {
        wu[i - 1] = wa[i - 1] * wu[i] + wb[i - 1];
        wv[i - 1] = wa[i - 1] * wv[i] + wc[i - 1];
    }

    FuncType k = (wb[N] + wa[N] * wu[1]) / (1 - wc[N] - wa[N] * wv[1]);
    FuncType c1tmp = k;
    for(size_t i = N; i > 0; --i)
    {
        c1[i - 1] = wu[i - 1] + k * wv[i - 1];
        // Optimization
        VarType h = x[i] - x[i - 1];
        FuncType v = f[i] - f[i - 1];
        // Calculate remaining coefficients
        c2[i - 1] = (3 * v - h * (2 * c1[i - 1] + c1tmp))
                     / (h * h);
        c3[i - 1] = (-2 * v + h * (c1[i - 1] + c1tmp))
                     / (h * h * h);
        c1tmp = c1[i - 1];
    }
}
