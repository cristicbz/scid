module scid.splines.univariate.cubic;

import scid.common.meta;
import scid.splines.univariate.base;

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
    /* TODO: After testing, reduce the number of workspaces by using
     *       the coefficient arrays as workspaces.
     *       Attention!
     *       If FuncType operations are considerably slower than VarType operations
     *       _wa, _wc and all _aX arryas should be of type VarType. So, no storage
     *       in coefficient arrays for them.
     */

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

        // Workspaces
        VarType _wa[];
        FuncType _wb[];
        VarType _wc[];
        FuncType _wu[];
        FuncType _wv[];

        void _allocContents(size_t maxSize) // TODO: use scid.core.memory
        {
            _c1.length = maxSize;
            _c2.length = maxSize;
            _c3.length = maxSize;
            _wa.length = maxSize;
            _wb.length = maxSize;
            _wc.length = maxSize;
            _wu.length = maxSize;
            _wv.length = maxSize;
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
            static assert(false, "No implementation for fixVar mode");

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
                // TODO: optimaize memory usage after testing of the algorithm
                if(_bcLeftType == BoundCond.periodic)
                {
                    // The index of the last point of spline
                    size_t N = pointsNum - 1;
                    // Some variables for optimization
                    VarType h0 = _x[1] - _x[0];
                    FuncType v0 = _f[1] - _f[0];
                    FuncType d0 = v0 / h0;
                    // Linear equation coefficients
                    VarType a0, a1, a2;
                    // Linear equation right part
                    FuncType b;

                    // Boundary conditions on the left side (the first equation)
                    // TODO: optimize
                    // FIXME: check
                    VarType hf = _x[pointsNum-1] - _x[pointsNum-2];
                    a0 = 1 / hf;
                    a2 = 1 / h0;
                    a1 = 2 * (a0 + a2);
                    b = 3 * ((_f[1] - _f[0]) / (h0 * h0)
                             + (_f[0] - _f[pointsNum-2]) / (hf * hf));

                    // Build'n'sweep forward

                    // First step
                    _wa[0] = 0;
                    _wb[0] = Zero!(FuncType);
                    _wc[0] = 1;
                    // TODO: _wa -> _c3, wb -> _c1 _wc -> c2
                    for(size_t i = 1; i < N; ++i)
                    {
                        // Optimization
                        VarType h1 = _x[i+1] - _x[i];
                        FuncType v1 = _f[i+1] - _f[i];
                        FuncType d1 = v1 / h1;
                        // Calculate eqaution coefficients
                        a0 = h1 / (h0 + h1);
                        a1 = 2;
                        a2 = 1 - a0;
                        b = 3 * (a0 * d0 + a2 * d1);
                        // Forward sweep step
                        VarType factor = 1 / (a1 + a0 * _wa[i - 1]);
                        _wa[i] = -a2 * factor;
                        _wb[i] = (b - a0 * _wb[i - 1]) * factor;
                        _wc[i] = -a0 * _wc[i - 1] * factor;
                        // Optimization
                        h0 = h1;
                        v0 = v1;
                        d0 = d1;
                    }

                    // Boundary conditions on the right side (the last equation)
                    a0 = 0;
                    a1 = 1;
                    a2 = -1;
                    b = Zero!(FuncType);

                    // TODO: optimize
                    VarType factor = 1 / (a1 + a0 * _wa[N - 1]);
                    _wa[N] = -a2 * factor;
                    _wb[N] = (b - a0 * _wb[N - 1]) * factor;
                    _wc[N] = -a0 * _wc[N - 1] * factor;

                    // Sweep backward
                    // TODO: _wu -> _c2, wv -> _c3
                    _wu[N - 1] = _wb[N-1];
                    _wv[N - 1] = One!(FuncType)*(_wa[N-1] + _wc[N-1]);
                    for(size_t i = N - 1; i > 0; --i)
                    {
                        _wu[i - 1] = _wa[i] * _wu[i] + _wb[i];
                        _wv[i - 1] = _wa[i] * _wv[i] + _wc[i];
                    }

                    FuncType k = (_wb[N] + _wa[N] * _wu[1])
                              / (1 - _wc[N] - _wa[N] * _wv[1]);
                    _c1[N] = k;
                    for(size_t i = N; i > 0; --i)
                    {
                        _c1[i - 1] = _wu[i - 1] + k * _wv[i - 1];
                        // Optimization
                        VarType h = _x[i] - _x[i - 1];
                        FuncType v = _f[i] - _f[i - 1];
                        // Calculate remaining coefficients
                        _c2[i - 1] = (3 * v - h * (2 * _c1[i - 1] + _c1[i]))
                                     / (h * h);
                        _c3[i - 1] = (-2 * v + h * (_c1[i - 1] + _c1[i]))
                                     / (h * h * h);
                    }
                }
                else
                {
                    // The index of the last point of spline
                    size_t N = pointsNum - 1;
                    // Some variables for optimization
                    VarType h0 = _x[1] - _x[0];
                    FuncType v0 = _f[1] - _f[0];
                    FuncType d0 = v0 / h0;
                    // Linear equation coefficients
                    VarType a0, a1, a2;
                    // Linear equation right part
                    FuncType b;

                    // Boundary conditions on the left side (the first equation)
                    if(_bcLeftType == BoundCond.deriv1)
                    {// FIXME: check
                        a1 = 1;
                        a2 = 0;
                        b = _bcLeftVal;
                    }
                    else if(_bcLeftType == BoundCond.deriv2)
                    {// FIXME: check
                        a1 = 2;
                        a2 = 1;
                        b = 3 * d0 - h0 * _bcLeftVal / 2;
                    }
                    // a0 won't be used in the first equation

                    // Build'n'sweep forward

                    // First step
                    VarType factor = 1 / a1;
                    _wa[0] = -a2 * factor;
                    _wb[0] = b * factor;

                    for(size_t i = 1; i < N; ++i)
                    {
                        // Optimization
                        VarType h1 = _x[i+1] - _x[i];
                        FuncType v1 = _f[i+1] - _f[i];
                        FuncType d1 = v1 / h1;
                        // Calculate eqaution coefficients
                        a0 = h1 / (h0 + h1);
                        a1 = 2;
                        a2 = 1 - a0;
                        b = 3 * (a0 * d0 + a2 * d1);
                        // Forward sweep step
                        factor = 1 / (a1 + a0 * _wa[i - 1]);
                        _wa[i] = -a2 * factor;
                        _wb[i] = (b - a0 * _wb[i - 1]) * factor;
                        // Optimization
                        h0 = h1;
                        v0 = v1;
                        d0 = d1;
                    }

                    // Boundary conditions on the right side (the last equation)
                    if(_bcRightType == BoundCond.deriv1)
                    {// FIXME: check
                        a0 = 0;
                        a1 = 1;
                        b = _bcRightVal;
                    }
                    else if(_bcRightType == BoundCond.deriv2)
                    {// FIXME: check
                        a0 = 1;
                        a1 = 2;
                        b = 3 * d0 - h0 * _bcRightVal / 2;
                    }

                    // Sweep backward and calculate the coefficients
                    factor = 1 / (a1 + a0 * _wa[N - 1]);
                    _c1[N] = (b - a0 * _wb[N - 1]) * factor;
                    for(size_t i = N; i > 0; --i)
                    {
                        // Backward sweep step
                        _c1[i - 1] = _wa[i - 1] * _c1[i] + _wb[i - 1];
                        // Optimization
                        VarType h = _x[i] - _x[i - 1];
                        FuncType v = _f[i] - _f[i - 1];
                        // Calculate remaining coefficients
                        _c2[i - 1] = (3 * v - h * (2 * _c1[i - 1] + _c1[i]))
                                     / (h * h);
                        _c3[i - 1] = (-2 * v + h * (_c1[i - 1] + _c1[i]))
                                     / (h * h * h);
                    }
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
                        MySpline.BoundCond.deriv1,
                        MySpline.BoundCond.deriv1,
                        0, 27);
    assert(spl._calcFunction(2, 1) == 8);
    assert(spl._calcDeriv(2, 1) == 12);
}
