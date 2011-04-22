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

/** Structure to wrap a one-dimensional spline.
  * It caches the previous value of an index lookup. When the subsequent
  * interpolation point falls in the same interval its index value can be
  * returned immediately.
  * Many Spline1dView structures can refer the same spline.
  * Params:
  *     Tspline = type of spline to wrap
  */
struct SplineView(Tspline)
{
private:
    Tspline* _spline;

    alias typeof(_spline._f[0]) Tfunc;
    alias typeof(_spline._x[0]) Tvar;

    size_t _index = 0; // The index of current interval

    // Find the value x in a sorted array a
    void _updateIndex(Tvar x)
    {
        // Does current interval contain x?
        if((_spline._x[_index + 1] < x) || (_spline._x[_index] > x))
            // Find new interval if it doesn't
            _index = binarySearch(_spline._x, x);
    }

public:
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

/** The type of optimization */
enum SplineOptim
{
    /// no special features
    normal,
    /// accelerates multiple calculations for the same variable value array
    fixVar,
}

/** The type of the storage of the variable and function value arrays */
enum SplineStorage
{
    /// Do not make any copies
    none,
    /// Keep a copy of the variable value array
    var,
    /// Keep a copy of both variable and function value arrays
    all
}

/** Linear one-dimensional spline (order = 1, defect = 1).
  * Params:
  *     Tfunc = type of function
  *     Tvar = type of variable
  *     storage = SplineStorage.none by default
  *     optim = SplineOptim.normal by default
  */
struct SplineLinear(Tvar, Tfunc,
                    SplineOptim optim = SplineOptim.normal,
                    SplineStorage storage = SplineStorage.none)
{
private:
    Tvar[] _x; // variable value array
    Tfunc[] _f; // function value array

    // Spline parameters:
    Tfunc[] _c1;
    /* The interpolant is:
     *     f(x) = _f[i] + _c1[i] * dx
     *     dx = x - _x[i]
     */

    static if(storage != SplineStorage.none)
    {
        /* NOTE: Sometimes it's necessary to keep more memory for the point
         *       arrays, e.g. to add points in interactive applications.
         *       Buffers and slices looks like a better idea than _size member
         *       especially since _x.length is widely used as current size
         *       of the spline and not only in this structure.
         *       The member _buff is added just for consistency.
         */
        Tvar[] _bufx; // memory reserved for variable value array
        static if(storage == SplineStorage.all)
        {
            Tfunc[] _buff; // memory reserved for function value array
        }
    }

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
    }

    /* This is the place for numerical algorithms */
    static if(optim == SplineOptim.fixVar)
    {
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

public:
    /// Minimal number of points needed for the spline
    enum size_t minPoints = 2;

    this(size_t maxSize)
    {
        _allocContents(maxSize);
    }

    this(Tvar[] x, Tfunc[] y, bool calcNow = true)
    {
        this(x.length);
        setAll(x, y);
        if(calcNow)
            calculate();
    }

    void setVar(Tvar[] x)
    {
        static if(storage != SplineStorage.none)
        {
            if(x.length > _bufx.length)
                _allocContents(x.length);
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

    void setFunc(Tfunc[] f)
    {
        static if(storage == SplineStorage.all)
        {
            if(f.length > _buff.length)
                _allocContents(f.length);
            _f = _buff[0..f.length];
            _f[] = f[];
        }
        else
        {
            _f = f;
        }
    }

    void setAll(Tvar[] x, Tfunc[] f)
    {
        setVar(x);
        setFunc(f);
    }

    void calculate()
    in
    {
        assert(_f.length == _x.length,
               "Variable and function value arrays have different sizes");
        assert(_x.length >= minPoints, "Not enough points for interpolation");
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

    void calculate(Tvar[] x, Tfunc[] y)
    {
        setAll(x, y);
        calculate();
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
private:
    Tvar[] _x; // variable value array
    Tfunc[] _f; // function value array

    static if(storage != SplineStorage.none)
    {
        Tvar[] _bufx; // _x[] is a slice of _bufx[]

        static if(storage == SplineStorage.all)
        {
            Tfunc[] _buff; // _f[] is a slice of _buff[]
        }
    }

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
            // FIXME
        }
    }

    /* This is the place for numerical algorithms */
    static if(optim == SplineOptim.fixVar)
    {
        // TODO: implement this area when the algorithm will be tested

        bool _needUpdateVar;

        /* Data depending only on variable values */
        // FIXME

        void _calcVarDependent()
        {
            // FIXME
        }

        void _calcFuncDependent()
        {
            // FIXME
        }
    }
    else
    {
        void _calcAll()
        {
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
                // TODO: use _cX arrays as workspaces

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

    // Calculate function in a given interval
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

public:
    /// Minimal number of points needed for the spline
    enum size_t minPoints = 3;

    public {/** Boundary conditions (BC) for cubic spline. */
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

    this(size_t maxSize)
    {
        _allocContents(maxSize);
    }

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

    void setVar(Tvar[] x)
    {
        static if(storage != SplineStorage.none)
        {
            if(x.length > _bufx.length)
                _allocContents(x.length);
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

    void setFunc(Tfunc[] f)
    {
        static if(storage == SplineStorage.all)
        {
            if(f.length > _buff.length)
                _allocContents(f.length);
            _f = _buff[0..f.length];
            _f[] = f[];
        }
        else
        {
            _f = f;
        }
    }

    void setAll(Tvar[] x, Tfunc[] f)
    {
        setVar(x);
        setFunc(f);
    }

    void calculate()
    in
    {
        assert(_f.length == _x.length,
               "Variable and function value arrays have different sizes");
        assert(_x.length >= minPoints, "Not enough points for interpolation");
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

    void calculate(Tvar[] x, Tfunc[] y)
    {
        setAll(x, y);
        calculate();
    }
}
// TODO: write spline building algorithms as separate (non-member) functions
// TODO: unittest

struct SplineCubicA(Tfunc, Tvar) // FIXME: to be removed
{
private:
    Tfunc[] _f; // function value array
    Tvar[] _x; // variable value array

    // Spline parameters:
    Tfunc[] _c1;
    Tfunc[] _c2;
    Tfunc[] _c3;
    /* The interpolant is:
     *     f(x) = _f[i] + _c1[i] * dx + _c2[i] * dx*dx + _c3[i] * dx*dx*dx
     *     dx = x - _x[i]
     */

    // Calculate function in a given interval
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

public:
    /// Flags used by spline creation function
    enum flags
    {
        dupSource = 0x00_00_00_01, // duplicate point arrays
        update = 0x00_00_00_02, // calculate spline coefficients immediately
        bcPeriodic = 0x00_00_00_04,    /* periodic spline boundary condition
                                        * parameters are ignored */
        bcLeftDeriv1 = 0x00_00_00_08,  /* 1st BC parameter is the
                                        * 1st derivative at x[0] */
        bcLeftDeriv2 = 0x00_00_00_10,  /* 1st BC parameter is the
                                        * 2nd derivative at x[0] */
        bcRightDeriv1 = 0x00_00_00_20, /* 2nd BC parameter is the
                                        * 1st derivative at x[N] */
        bcRightDeriv2 = 0x00_00_00_40  /* 2nd BC parameter is the
                                        * 2nd derivative at x[N] */
    }
    enum flags defaultFlags = flags.update | flags.bcLeftDeriv2
                              | flags.bcRightDeriv2;
    /// Minimal number of points needed for the spline
    enum size_t minPoints = 3;

    this(Tfunc[] f, Tvar[] x, flags flg = defaultFlags,
         Tfunc bcLeft = 0, Tfunc bcRight = 0)
    {
        bind(f, x, flg, bcLeft, bcRight);
    }

    /** Select point arrays
      * Params:
      *     f = array of function values
      *     x = array of variable values
      *         must be sorted in ascending order if update flag is on
      */
    void bind(Tfunc[] f, Tvar[] x, flags flg = defaultFlags,
                        Tfunc bcLeft = 0, Tfunc bcRight = 0)
    in
    {
        assert(f.length == x.length,
               "Variable and function value arrays have different sizes");
        assert(x.length >= minPoints, "Not enough points for interpolation");
        for(size_t i = 0; i < x.length-1; ++i)
            assert(x[i] < x[i+1],
                   "Variable value array is not sorted in ascending order");
    }
    body
    {
        if(flg & flags.dupSource) // Make a copy of source data if required
        {
            _f = f.dup;
            _x = x.dup;
        }
        else // Otherwise just save links and hope it won't change
        {
            /* NOTE: will the garbage collector free old arrays now
                     if they were duplicates? */
            _f = f;
            _x = x;
        }
        if(flg & flags.update)
            update(flg, bcLeft, bcRight);
    }

    /** Calculate spline prarmeters
      */
    void update(flags flg = defaultFlags, Tfunc bcLeft = 0, Tfunc bcRight = 0)
    in
    {
        // FIXME: what if arrays are not allocated yet?
        // TODO: mixins for repeating code
        assert(_f.length == _x.length,
               "Variable and function value arrays have different sizes");
        assert(_x.length >= minPoints, "Not enough points for interpolation");
        for(size_t i = 0; i < _x.length-1; ++i)
            assert(_x[i] < _x[i+1],
                   "Variable value array is not sorted in ascending order");
    }
    body // FIXME: it needs more comments
    {
        // Tridiagonal matrix
        Tfunc[][3] a = [new Tfunc[_x.length],
                        new Tfunc[_x.length],
                        new Tfunc[_x.length]];
        // Right part of equation
        Tfunc[] b = new Tfunc[_x.length];

        if(flg & flags.bcPeriodic)
        {
            Tvar h0 = _x[1] - _x[2];
            Tvar hf = _x[$-1] - _x[$-2];
            a[0][0] = 1 / hf;
            a[1][0] = 2 * (1 / h0 + 1 / hf);
            a[2][0] = 1 / h0;
            b[0] = 3 * ((_f[1] - _f[0]) / (h0 * h0)
                        + (_f[0] - _f[$-2]) / (hf * hf));
        }
        else if(flg & flags.bcLeftDeriv1)
        {
            a[1][0] = 1;
            a[2][0] = 0;
            b[0] = bcLeft;
        }
        else if(flg & flags.bcLeftDeriv2)
        {
            a[1][0] = 2;
            a[2][0] = 1;
            Tvar h0 = _x[1] - _x[0];
            b[0] = 3 * (_f[1] - _f[0]) / h0 - h0 * bcLeft / 2;
        }

        Tvar h0 = _x[1] - _x[0];
        Tfunc v0 = _f[1] - _f[0];
        Tfunc d0 = v0 / h0;
        for(size_t i = 1; i < _x.length - 1; ++i)
        {
            Tvar h1 = _x[i+1] - _x[i];
            Tfunc v1 = _f[i+1] - _f[i];
            Tfunc d1 = v1 / h1;

            a[0][i] = h0 / (h0 + h1);
            a[1][i] = 2;
            a[2][i] = 1 - a[0][i];
            b[i] = 3 * (a[0][i] * d1 + a[2][i] * d0);

            h0 = h1;
            v0 = v1;
            d0 = d1;
        }

        if(flg & flags.bcPeriodic)
        {
            a[0][$-1] = 0;
            a[1][$-1] = 1;
            a[2][$-1] = -1;
            b[$-1] = 0;
        }
        else if(flg & flags.bcRightDeriv1)
        {
            a[0][$-1] = 0;
            a[1][$-1] = 1;
            b[$-1] = bcRight;
        }
        else if(flg & flags.bcRightDeriv2)
        {
            a[0][$-1] = 1;
            a[1][$-1] = 2;
            Tvar hf = _x[$-1] - _x[$-2];
            b[$-1] = 3 * (_f[$-1] - _f[$-2]) / hf - hf * bcRight / 2;
        }

        _c1.length = _x.length;
        if(flg & flags.bcPeriodic)
            SweepCycl3Diag(a, b, _c1);
        else
            Sweep3Diag(a, b, _c1);

        _c2.length = _x.length - 1;
        _c3.length = _x.length - 1;

        for(size_t j = _x.length - 1; j > 0; --j)
        {
            size_t i = j - 1; // NOTE: because there is no -1 for size_t
            Tvar h = _x[i+1] - _x[i];
            Tfunc v = _f[i+1] - _f[i];
            _c2[i] = (3 * v - h * (2 * _c1[i] + _c1[i+1])) / (h * h);
            _c3[i] = (-2 * v + h * (_c1[i] + _c1[i+1])) / (h * h * h);
        }
    }
}

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

/* Solve Ax = b linear system for three-diagonal matrix A
 * a[1][0] * x[0] + a[2][0] * x[1] = b[0]
 * a[0][i] * x[i - 1] + a[1][i] * x[i] + a[2][i] * x[i + 1] = b[i],
 * i = 1 .. (N - 1)
 * a[0][N] * x[N - 1] + a[1][N] * x[N] = b[N] */
/* NOTE: not tested properly */
void sweep3Diag(T)(in T[][3] a, in T[] b, T[] x)
in
{ // FIXME: make assert strings more    informative
    assert(b.length >= 3, "b");
    assert(a[0].length == b.length, "a0");
    assert(a[1].length == b.length, "a1");
    assert(a[2].length == b.length, "a2");
    assert(x.length == b.length, "x");
}
body
{
    size_t N = b.length - 1; // Index of the last element in x array
    /* FIXME: the workspaces should be passed as parameters
              or allocated using scid.core.memory */
    T[] wsa = new T[N + 1];
    T[] wsb = new T[N + 1];

    // forward sweep
    wsa[0] = -a[2][0] / a[1][0];
    wsb[0] = b[0] / a[1][0];
    for(size_t i = 1; i <= N; ++i)
    {
        T factor = 1 / (a[1][i] + a[0][i] * wsa[i - 1]);
        wsa[i] = -a[2][i] * factor; // NOTE: wsa[N] is never used
        wsb[i] = (b[i] - a[0][i] * wsb[i - 1]) * factor;
    }

    // backward sweep
    x[N] = wsb[N];
    for(size_t i = N; i > 0; --i)
        x[i - 1] = wsa[i - 1] * x[i] + wsb[i - 1];
}

/* Solve Ax = b linear system for almost three-diagonal matrix A
 * using cyclic sweep
 * a[0][0] * x[N - 1] + a[1][0] * x[0] + a[2][0] * x[1] = b[0]
 * a[0][i] * x[i - 1] + a[1][i] * x[i] + a[2][i] * x[i + 1] = b[i],
 * i = 1 .. (N - 1) */
/* NOTE: not tested properly */
void sweepCycl3Diag(T)(in T[][3] a, in T[] b, T[] x)
in
{ // FIXME: make assert strings more    informative
    assert(b.length >= 3, "b");
    assert(a[0].length == b.length, "a0");
    assert(a[1].length == b.length, "a1");
    assert(a[2].length == b.length, "a2");
    assert(x.length == b.length, "x");
}
body
{
    size_t N = b.length - 1; // Index of the last element in x array
    /* FIXME: the workspaces should be passed as parameters
              or allocated using scid.core.memory */
    T[] wsa = new T[N + 1];
    T[] wsb = new T[N + 1];
    T[] wsc = new T[N + 1];
    T[] wsv = new T[N + 1];

    /* forward sweep
     * wsa[i - 1] is alpha[i]
     * wsb[i - 1] is beta[i]
     * wsc[i - 1] is gamma[i] */
    wsa[0] = 0;
    wsb[0] = 0;
    wsc[0] = 1;
    for(size_t i = 1; i <= N; ++i)
    {
        T factor = 1 / (a[1][i] + a[0][i] * wsa[i - 1]);
        wsa[i] = -a[2][i] * factor;
        wsb[i] = (b[i] - a[0][i] * wsb[i - 1]) * factor;
        wsc[i] = -a[0][i] * wsc[i - 1] * factor;
    }

    /* forward sweep
     * x[i] is u[i]
     * wsv[i] is v[i]
     */
    x[N - 1] = wsb[N - 1];
    wsv[N - 1] = wsa[N - 1] + wsc[N - 1];
    for(size_t i = N - 2; i > 0; --i)
    {
        x[i] = wsa[i + 1] * x[i + 1] + wsb[i + 1];
        wsv[i] = wsa[i + 1] * wsv[i + 1] + wsc[i + 1];
    }

    /* final step
     * x[i] is y[i]
     */
    x[0] = (wsb[N] + wsa[N] * x[1]) / (1 - wsc[N] - wsa[N] * wsv[1]);
    for(size_t i = 1; i < N; ++i )
        x[i] += x[0] * wsv[i];
    x[N] = x[0];
}

