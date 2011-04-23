/** This module contains types for interpolation.
  *
  * This is an alpha version. The interface can be changed.
  *
  * Version: 0.2-a
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev. All rights reserved.
  * License: Boost License 1.0
  */
module scid.interpolation;
// TODO: add more embedded documentation
// FIXME: add "in" and "out" qualifiers for function arguments

/** Determines which arras should be duplicated and stored in the spline
  * structue.
  */
enum SplineStorage
{
    none, /// Do not make any copies
    var, /// Keep a copy of the variable value array
    all /// Keep a copy of both variable and function value arrays
}

/// The type of optimization
enum SplineOptim
{
    // No, no, no! David Blaine! We don't need your magic!
    normal, /// no special features
    // Use the Force, Luke!
    fixVar, /** accelerate multiple calculations with the same
              * variable value array
              */
}

/** Common functions for univariate splines
  */
mixin template splineBase(Tvar, Tfunc,
                          SplineStorage storage)
{
    // Variable and function value arrays
    private
    {
        size_t _maxSize = 0;
        Tvar[] _x; // variable value array
        Tfunc[] _f; // function value array

        static if(storage != SplineStorage.none)
        {
            Tvar[] _bufx; // memory reserved for variable value array
            static if(storage == SplineStorage.all)
            {
                Tfunc[] _buff; // memory reserved for function value array
            }
        }
    }

    public
    {
        /** Select a variable value array for the spline
          */
        void setVar(Tvar[] x)
        in
        {
            static if(storage != SplineStorage.none)
            {
                assert(x.length <= _maxSize,
                       "Variable value array is too large");
            }
        }
        body
        {
            static if(storage != SplineStorage.none)
            {
                _x = _bufx[0..x.length];
                _x[] = x[];
            }
            else
            {
                _x = x;
            }
            static if(optim == SplineOptim.fixVar)
            {
                _needUpdateVar = true;
            }
        }

        /** Select a function value array for the spline
          */
        void setFunc(Tfunc[] f)
        in
        {
            static if(storage == SplineStorage.all)
            {
                assert(f.length <= _maxSize,
                       "Function value array is too large");
            }
        }
        body
        {
            static if(storage == SplineStorage.all)
            {
                _f = _buff[0..f.length];
                _f[] = f[];
            }
            else
            {
                _f = f;
            }
        }

        /** Select variable and function value arrays for the spline
          */
        void setAll(Tvar[] x, Tfunc[] f)
        {
            // TODO: enlarge buffers if they are too small
            setVar(x);
            setFunc(f);
        }

        /** Calculate the interpolant data. Makes it ready for function evaluation.
          */
        void calculate()
        in
        {
            assert(_f.length == _x.length,
                   "Variable and function value arrays have different sizes");
            assert(_x.length >= minPoints,
                   "Not enough points for interpolation");
            for(size_t i = 0; i < _x.length - 1; ++i)
                assert(_x[i] < _x[i+1],
                       "Variable value array is not sorted in ascending order");
        }
        body
        {
            static if(optim == SplineOptim.fixVar)
            {
                if(_needUpdateVar)
                {
                    _calcVarDependent();
                    _needUpdateVar = false;
                }
                _calcFuncDependent();
            }
            else
            {
                _calcAll();
            }
        }

        /** Calculate the interpolant data for given variable and function value
          * arrays. Makes it ready for function evaluation.
          */
        void calculate(Tvar[] x, Tfunc[] y)
        {
            setAll(x, y);
            calculate();
        }
    }
}

/** Structure to access a one-dimensional spline, which provides the interface
  * for function evaluation and accelerates this action if multiple successive
  * evaluations are made on the same interval.
  *
  * Many SplineView structures can refer the same spline.
  *
  * Under the hood, it caches the previous value of an index lookup and return
  * it immediately if the subsequent interpolation point falls in the same
  * interval.
  *
  * Params:
  *     Tspline = type of spline to wrap
  *
  * Examples:
  * --------------
  * double[] x = [0.0, 1.0];
  * double[] y = [1.0, 3.0];
  * auto spline = SplineLinear!(double, double)(x, y);
  * auto splineView = SplineView!(typeof(spline))(spline);
  * dobule f = splineView.eval(0.4);
  * --------------
  */
struct SplineView(Tspline)
{
    // Pointer to the wrapped spline
    private Tspline* _spline; // TODO: make it a property?

    // Types from the wrapped spline
    public
    {
        /** Function and variable types for the wrapped spline */
        alias typeof(_spline._f[0]) Tfunc;
        alias typeof(_spline._x[0]) Tvar; ///ditto
    }

    // Index lookup
    private
    {
        size_t _index = 0; // The index of current interval

        // Find the value x in a sorted array a
        void _updateIndex(Tvar x)
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
        Tfunc eval(Tvar x)
        in
        {
            assert((x >= _spline._x[0]) && (x <= _spline._x[$-1]),
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
        Tfunc deriv(Tvar x)
        in
        {
            assert((x >= _spline._x[0]) && (x <= _spline._x[$-1]),
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

/** Linear one-dimensional spline (order = 1, defect = 1).
  *
  * Params:
  *     Tfunc = type of function
  *     Tvar = type of variable
  *     storage = spline storage type (none by default)
  *     optim = spline optimization type (normal by default)
  */
struct SplineLinear(Tvar, Tfunc,
                    SplineOptim optim = SplineOptim.normal,
                    SplineStorage storage = SplineStorage.none)
{
    mixin splineBase!(Tvar, Tfunc, storage);

    // Data
    private
    {
        // Spline parameters:
        Tfunc[] _c1;
        /* The interpolant is:
         *     f(x) = _f[i] + _c1[i] * dx
         *     dx = x - _x[i]
         */

        void _allocContents(size_t maxSize) // TODO: use scid.core.memory
        {
            static if(storage != SplineStorage.none)
            {
                _bufx.length = maxSize;
                static if(storage == SplineStorage.all)
                {
                    _buff.length = maxSize;
                }
            }
            _c1.length = maxSize - 1;
            static if(optim == SplineOptim.fixVar)
            {
                _rdx.length = maxSize - 1;
            }
            _maxSize = maxSize;
        }
    }

    // Numerical core
    private
    {
        static if(optim == SplineOptim.fixVar)
        {
            /* NOTE: division is slower than multiplication.
             */
            bool _needUpdateVar;

            /* Data depending only on variable values */
            Tvar[] _rdx;

            void _calcVarDependent()
            {
                for(size_t i = 0; i < _x.length-1; ++i)
                    _rdx[i] = 1 / (_x[i+1] - _x[i]);
            }

            void _calcFuncDependent()
            {
                for(size_t i = 0; i < _x.length - 1; ++i)
                    _c1[i] = (_f[i+1] - _f[i]) * _rdx[i];
            }
        }
        else
        {
            void _calcAll()
            {
                for(size_t i = 0; i < _x.length - 1; ++i)
                    _c1[i] = (_f[i+1] - _f[i]) / (_x[i+1] - _x[i]);
            }
        }
    }

    // Function evaluation code
    private
    {
        // Calculate function in a given interval
        Tfunc _calcFunction(Tvar x, size_t index)
        {
            return _f[index] + _c1[index] * (x - _x[index]);
        }

        // Calculate first derivative in a given interval
        Tfunc _calcDeriv(Tvar x, size_t index)
        {
            return _c1[index];
        }
    }

    public
    {
        /// Minimal number of points needed for the spline
        enum size_t minPoints = 2;

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
        this(Tvar[] x, Tfunc[] y, bool calcNow = true)
        {
            this(x.length);
            setAll(x, y);
            if(calcNow)
                calculate();
        }
    }
} // TODO: unittest

/* -------------------------------------------------------------------------- */

/** Cubic one-dimensional spline (order = 3, defect = 1).
  * Params:
  *     Tfunc = type of function
  *     Tvar = type of variable
  */
struct SplineCubic(Tvar, Tfunc,
                   SplineOptim optim = SplineOptim.normal,
                   SplineStorage storage = SplineStorage.none)
{
    /* TODO: After testing, reduce the number of workspaces by using
     *       the coefficient arrays as workspaces.
     *       Attention!
     *       If Tfunc operations are considerably slower than Tvar operations
     *       _wa, _wc and all _aX arryas should be of type Tvar. So, no storage
     *       in coefficient arrays for them.
     */

    mixin splineBase!(Tvar, Tfunc, storage);

    // Data and workspaces
    private
    {
        // Spline parameters:
        Tfunc[] _c1;
        Tfunc[] _c2;
        Tfunc[] _c3;
        /* The interpolant is:
         *     f(x) = _f[i] + _c1[i] * dx + _c2[i] * dx*dx + _c3[i] * dx*dx*dx
         *     dx = x - _x[i]
         */

        // Workspaces
        Tvar _wa[];
        Tfunc _wb[];
        Tvar _wc[];
        Tfunc _wu[];
        Tfunc _wv[];

        void _allocContents(size_t maxSize) // TODO: use scid.core.memory
        {
            static if(storage != SplineStorage.none)
            {
                _bufx.length = maxSize;
                static if(storage == SplineStorage.all)
                {
                    _buff.length = maxSize;
                }
            }
            _c1.length = maxSize - 1;
            _c2.length = maxSize - 1;
            _c3.length = maxSize - 1;
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
    private
    {
        static if(optim == SplineOptim.fixVar)
        {
            // TODO: implement this area when the algorithm will be tested

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
                    size_t N = _x.length - 1;
                    // Some variables for optimization
                    Tvar h0 = _x[1] - _x[0];
                    Tfunc v0 = _f[1] - _f[0];
                    Tfunc d0 = v0 / h0;
                    // Linear equation coefficients
                    Tvar a0, a1, a2;
                    // Linear equation right part
                    Tfunc b;

                    // Boundary conditions on the left side (the first equation)
                    // TODO: optimize
                    Tvar hf = _x[$-1] - _x[$-2];
                    a0 = 1 / hf;
                    a2 = 1 / h0;
                    a1 = 2 * (a0 + a2);
                    b = 3 * ((_f[1] - _f[0]) / (h0 * h0)
                             + (_f[0] - _f[$-2]) / (hf * hf));

                    // Build'n'sweep forward

                    // First step
                    _wa[0] = 0;
                    _wb[0] = 0;
                    _wc[0] = 1;
                    // TODO: _wa -> _c3, wb -> _c1 _wc -> c2
                    for(size_t i = 1; i < N; ++i)
                    {
                        // Optimization
                        Tvar h1 = _x[i+1] - _x[i];
                        Tfunc v1 = _f[i+1] - _f[i];
                        Tfunc d1 = v1 / h1;
                        // Calculate eqaution coefficients
                        a0 = h0 / (h0 + h1);
                        a1 = 2;
                        a2 = 1 - a0;
                        b = 3 * (a0 * d1 + a2 * d0);
                        // Forward sweep step
                        Tvar factor = 1 / (a1 + a0 * _wa[i - 1]);
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
                    b = 0;

                    // TODO: optimize
                    Tvar factor = 1 / (a1 + a0 * _wa[N - 1]);
                    _wa[N] = -a2 * factor;
                    _wb[N] = (b - a0 * _wb[N - 1]) * factor;
                    _wc[N] = -a0 * _wc[N - 1] * factor;

                    // Sweep backward
                    // TODO: _wu -> _c2, wv -> _c3
                    _wu[N - 1] = _wb[N-1];
                    _wv[N - 1] = _wa[N-1] + _wc[N-1];
                    for(size_t i = N - 1; i > 0; --i)
                    {
                        _wu[i - 1] = _wa[i] * _wu[i] + _wb[i];
                        _wv[i - 1] = _wa[i] * _wv[i] + _wc[i];
                    }

                    Tfunc k = (_wb[N] + _wa[N] * _wu[1])
                              / (1 - _wc[N] - _wa[N] * _wv[1]);
                    _c1[N] = k;
                    for(size_t i = N; i > 0; --i)
                    {
                        _c1[i - 1] = _wu[i - 1] + k * _wv[i - 1];
                        // Optimization
                        Tvar h = _x[i] - _x[i - 1];
                        Tfunc v = _f[i] - _f[i - 1];
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
                    size_t N = _x.length - 1;
                    // Some variables for optimization
                    Tvar h0 = _x[1] - _x[0];
                    Tfunc v0 = _f[1] - _f[0];
                    Tfunc d0 = v0 / h0;
                    // Linear equation coefficients
                    Tvar a0, a1, a2;
                    // Linear equation right part
                    Tfunc b;

                    // Boundary conditions on the left side (the first equation)
                    if(_bcLeftType == BoundCond.deriv1)
                    {
                        a1 = 1;
                        a2 = 0;
                        b = _bcLeftVal;
                    }
                    else if(_bcLeftType == BoundCond.deriv2)
                    {
                        a1 = 2;
                        a2 = 1;
                        b = 3 * d0 - h0 * _bcLeftVal / 2;
                    }
                    // a0 won't be used in the first equation

                    // Build'n'sweep forward

                    // First step
                    Tfunc factor = 1/ a1;
                    _wa[0] = -a2 * factor;
                    _wb[0] = b * factor;

                    for(size_t i = 1; i < N; ++i)
                    {
                        // Optimization
                        Tvar h1 = _x[i+1] - _x[i];
                        Tfunc v1 = _f[i+1] - _f[i];
                        Tfunc d1 = v1 / h1;
                        // Calculate eqaution coefficients
                        a0 = h0 / (h0 + h1);
                        a1 = 2;
                        a2 = 1 - a0;
                        b = 3 * (a0 * d1 + a2 * d0);
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
                    if(_bcLeftType == BoundCond.deriv1)
                    {
                        a0 = 0;
                        a1 = 1;
                        b = _bcRightVal;
                    }
                    else if(_bcLeftType == BoundCond.deriv2)
                    {
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
                        Tvar h = _x[i] - _x[i - 1];
                        Tfunc v = _f[i] - _f[i - 1];
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
    private
    {
        Tfunc _calcFunction(Tvar x, size_t index)
        {
            double dx = x - _x[index];
            return _f[index] + dx * (_c1[index]
                                     + dx * (_c2[index] + dx * _c3[index]));
        }

        // Calculate first derivative in a given interval
        Tfunc _calcDeriv(Tvar x, size_t index)
        {
            double dx = x - _x[index];
            return _c1[index] + dx * (2 * _c2[index] + dx * 3 * _c3[index]);
        }
    }

    // Boundary conditions
    public
    {
        /// BC types
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

        BoundCond bcLeftType()
        {
            return _bcLeftType;
        }

        private Tfunc _bcLeftVal = 0;

        void bcLeftVal(Tfunc val)
        {
            _bcLeftVal = val;
        }

        Tfunc bcLeftVal()
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

        BoundCond bcRightType()
        {
            return _bcRightType;
        }

        private Tfunc _bcRightVal = 0;

        void bcRightVal(Tfunc val)
        {
            _bcRightVal = val;
        }

        Tfunc bcRightVal()
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
        this(Tvar[] x, Tfunc[] y, bool calcNow = true,
             BoundCond bcLType = BoundCond.deriv2,
             BoundCond bcRType = BoundCond.deriv2,
             Tfunc bcLVal = 0,
             Tfunc bcRVal = 0)
        {
            this(x.length);
            setAll(x, y);
            /* For the following case: this(x, y, calcNow, periodic) */
            if(bcLType == BoundCond.periodic)
                bcRType = BoundCond.periodic;
            if(calcNow)
                calculate();
        }
    }
} // TODO: unittest

/* -------------------------------------------------------------------------- */

// TODO: 1d Akima spline
// TODO: 1d B-spline
// TODO: non-spline interpolators
// TODO: 2d bilinear interpolation
// TODO: 2d bicubic interpolation

/* ========================================================================== */
/* NOTE: Perhaps these functions should go to a separte module */
private:

// Common binary search in a sorted array
// FIXME: Does phobos already have such function?
private pure size_t binarySearch(T)(T[] a, T x)
{
    size_t ilo = 0;
    size_t ihi = a.length;
    while(ihi - ilo > 1)
    {
        size_t i = (ihi + ilo) / 2;
        if(a[i] > x)
            ihi = i;
        else
            ilo = i;
    }
    if(ilo == a.length - 1)
        --ilo;
    return ilo;
}

