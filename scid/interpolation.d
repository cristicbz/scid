/** This module contains types for interpolation.
  *
  * VVA is variable value array
  * FVA is function value array
  *
  * This is an alpha version. The interface can be changed.
  *
  * Version: 0.3-a
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev. All rights reserved.
  * License: Boost License 1.0
  *
  */
module scid.interpolation;
// TODO: add more embedded documentation
// FIXME: add "in" and "out" qualifiers for function arguments

import std.math;
import std.algorithm;
import scid.common.traits;
import scid.vector;

/// The type of optimization
enum SplineOptim
{
    normal, /// No special features.
    fixVar /** Preliminary calculations depending only on VVA are made if it
             * changes. Minimal amount of operations is performed when FVA
             * changes.
             */
}

/** Common functions for univariate splines.
  * "Eoc" prefix means Element Or Container
  */
mixin template splineBase(EocVar, EocFunc)
{
    // Test the environment this template is mixed in
    static assert(is(typeof(optim) == SplineOptim),
                  "Spline must have \"optim\" parameter");
    static assert(is(typeof(this.minPoints) == size_t),
                  "Spline must have \"size_t minPoints\" field.");
    // NOTE: will not work because of a compiler bug
    // TODO: send a bug request
    /*
    static if(optim == SplineOptim.fixVar)
    {
        static assert(is(typeof(this._calcVarDependent) ==
                         void function(void)),
                      "Spline must have \"_calcVarDependent\" method");
        static assert(is(typeof(this._calcFuncDependent) ==
                         void function(void)),
                      "Spline must have \"_calcFuncDependent\" method");
    }
    else
    {
        static assert(is(typeof(this._calcAll) == void function(void)),
                      "Spline must have \"_calcAll\" method ");
    }*/

    alias ProduceArray!EocVar VarArray;
    alias ProduceArray!EocFunc FuncArray;
    alias BaseElementType!VarArray VarType;
    alias BaseElementType!FuncArray FuncType;

    private
    {
        size_t _maxSize = 0;
        VarArray _x; // variable value array (VVA)
        FuncArray _f; // function value array (FVA)

        static if(optim == SplineOptim.fixVar)
        {
            /* Flag indicating that all variable dependent data must be
               recalculated */
            bool _needUpdateVar;
        }
    }

    public
    {
        @property const length()
        {
            return _x.length;
        }

        /** Select a VVA for the spline
          */
        void setVar(VarArray x)
        {
            _x = x; // Copy source vector
            static if(optim == SplineOptim.fixVar)
            {
                _needUpdateVar = true; // Set flag since VVA has been changed
            }
        }

        /** Select a FVA for the spline
          */
        void setFunc(FuncArray f)
        {
            _f = f;
        }

        /** Select both VVA and FVA for the spline
          */
        void setAll(VarArray x, FuncArray f)
        {
            setVar(x);
            setFunc(f);
        }

        /** Calculate the interpolant data. Makes it ready for function
          * evaluation.
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
            /* FIXME: use isSorted when random-access interface
               will be implemented for vectors */
            /*assert(isSorted(_x),
                   "Variable value array is not sorted in ascending order");*/
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

        /** Calculate the interpolant data for given VVA and FVA. Makes it ready
          * for function evaluation.
          */
        void calculate(VarArray x, FuncArray y)
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
            // TODO: change _spline.length to $ when Vector implement opDollar
            assert((x >= _spline._x[0]) && (x <= _spline._x[_spline.length-1]),
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
            // TODO: change _spline.length to $ when Vector implement opDollar
            assert((x >= _spline._x[0]) && (x <= _spline._x[_spline.length-1]),
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
  *     EocVar = type of variable
  *     EocFunc = type of function
  *     optim = spline optimization type (normal by default)
  */
struct SplineLinear(EocVar, EocFunc,
                    SplineOptim optim = SplineOptim.normal)
{
    mixin splineBase!(EocVar, EocFunc);

    // Data
    private
    {
        // Spline parameters:
        FuncType[] _c1;
        /* The interpolant is:
         *     f(x) = _f[i] + _c1[i] * dx
         *     dx = x - _x[i]
         */

        void _allocContents(size_t maxSize) // TODO: use scid.core.memory
        {
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

            /* Data depending only on variable values */
            VarType[] _rdx;

            void _calcVarDependent()
            {
                for(size_t i = 0; i < _x.length - 1; ++i)
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
        FuncType _calcFunction(VarType x, size_t index)
        {
            return _f[index] + _c1[index] * (x - _x[index]);
        }

        // Calculate first derivative in a given interval
        FuncType _calcDeriv(VarType x, size_t index)
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
        this(VarArray x, FuncArray y, bool calcNow = true)
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
    private
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
                    size_t N = _x.length - 1;
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
                    VarType hf = _x[length-1] - _x[length-2];
                    a0 = 1 / hf;
                    a2 = 1 / h0;
                    a1 = 2 * (a0 + a2);
                    b = 3 * ((_f[1] - _f[0]) / (h0 * h0)
                             + (_f[0] - _f[length-2]) / (hf * hf));

                    // Build'n'sweep forward

                    // First step
                    _wa[0] = 0;
                    _wb[0] = 0;
                    _wc[0] = 1;
                    // TODO: _wa -> _c3, wb -> _c1 _wc -> c2
                    for(size_t i = 1; i < N; ++i)
                    {
                        // Optimization
                        VarType h1 = _x[i+1] - _x[i];
                        FuncType v1 = _f[i+1] - _f[i];
                        FuncType d1 = v1 / h1;
                        // Calculate eqaution coefficients
                        a0 = h0 / (h0 + h1);
                        a1 = 2;
                        a2 = 1 - a0;
                        b = 3 * (a0 * d1 + a2 * d0);
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
                    b = 0;

                    // TODO: optimize
                    VarType factor = 1 / (a1 + a0 * _wa[N - 1]);
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
                    size_t N = _x.length - 1;
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
                    FuncType factor = 1/ a1;
                    _wa[0] = -a2 * factor;
                    _wb[0] = b * factor;

                    for(size_t i = 1; i < N; ++i)
                    {
                        // Optimization
                        VarType h1 = _x[i+1] - _x[i];
                        FuncType v1 = _f[i+1] - _f[i];
                        FuncType d1 = v1 / h1;
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
    private
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

        private FuncType _bcLeftVal = 0;

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

        private FuncType _bcRightVal = 0;

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
        this(FuncArray x, VarArray y, bool calcNow = true,
             BoundCond bcLType = BoundCond.deriv2,
             BoundCond bcRType = BoundCond.deriv2,
             FuncType bcLVal = 0,
             FuncType bcRVal = 0)
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

/** One-dimensional Akima interpolation.
  *
  * Params:
  *     EocVar = type of variable
  *     EocFunc = type of function
  *     optim = spline optimization type (normal by default)
  *     Props = strings, describing properties of FuncType type to process.
  *             It is necessary because for this kind of splines not only
  *             linear operations are performed with function values.
  *
  * Examples:
  * Complex function
  * ----------
  * alias SplineAkima!(double, Complex!(double),
  *                    SplineOptim.normal,
  *                    ".re", ".im") MySpline;
  * ----------
  */
struct SplineAkima(EocVar, EocFunc,
                   SplineOptim optim = SplineOptim.normal,
                   Props...)
{
    // TODO: Add different boundary conditions
    // TODO: Add a mechanism for adding points to the curve
    // TODO: Implement support of compound types

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
    private
    {
        // Some templates for compound types
        private
        {
            string codeWeight(string w, string dl, string dr, Props...)()
            {
                string result = "";
                foreach(p; Props)
                {
                    static assert(is(typeof(p) == string),
                                  "Properties must be strings");
                    result ~= w ~ p ~
                              " = abs(" ~
                                  dl ~ p
                                  ~ " - " ~
                                  dr ~ p
                              ~ ");\n";
                }
                return result;
            }

            string codeCoeff1(string t, string wl, string wr,
                              string dl, string dr,
                              Props...)()
            {
                string result = "";
                foreach(p; Props)
                {
                    static assert(is(typeof(p) == string),
                                  "Properties should be strings");
                    /* Some template magic:
                     * ----------
                     * if((wl.p + wr.p) > 0)
                     *     t.p = (wl.p * d.p + wr.p * d.p) / (wl.p + wr.p);
                     * else
                     *     t.p = (dl.p + d.p) / 2;
                     * ----------
                     */
                    result ~= "if(" ~
                                  "(" ~
                                      wl ~ p
                                      ~ " + " ~
                                      wr ~ p
                                  ~ ")"
                                  ~ " > 0"
                              ~ ")\n    " ~
                                  t ~ p ~ " = " ~
                                  "(" ~
                                      wl ~ p
                                      ~ " * " ~
                                      dl ~ p
                                      ~
                                      " + "
                                      ~
                                      wr ~ p
                                      ~ " * " ~
                                      dr ~ p
                                  ~ ")"
                                  ~ " / " ~
                                  "(" ~
                                      wl ~ p
                                      ~ " + " ~
                                      wr ~ p
                                  ~ ");\n"
                              ~ "else\n    " ~
                                  t ~ p ~ " = " ~
                                  "(" ~
                                      dl ~ p
                                      ~ " + " ~
                                      dr ~ p
                                  ~ ")" ~
                                  " / 2;\n";
                }
                return result;
            }
        }

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
                /* TODO: After testing, refuse the workspaces by using
                 *       the coefficient arrays instead.
                 */
                size_t N = _x.length - 1;
                FuncType[] d = new FuncType[N + 2];
                FuncType[] w = new FuncType[N + 3];
                // Calculate slopes and weights
                for(size_t i = 1; i <= N; ++i)
                    d[i] = (_f[i] - _f[i - 1]) / (_x[i] - _x[i - 1]);
                d[0] = 2 * d[1] - d[2];
                d[N + 1] = 2 * d[N] - d[N - 1];
                for(size_t i = 1; i <= N + 1; ++i)
                    w[i] = abs(d[i] - d[i - 1]); // FIXME: real only
                w[0] = w[1];
                w[N + 2] = w[N + 1];
                for(size_t i = 0; i <= N; ++i)
                    if(w[i] + w[i + 2] > 0)
                        _c1[i] = (w[i] * d[i] + w[i + 2] * d[i + 1])
                                 / (w[i] + w[i + 2]);
                    else
                        _c1[i] = (d[i] + d[i + 2]) / 2;
                // Calculate the remaining coefficients
                for(size_t i = 0; i < N; ++i)
                {
                    VarType h = _x[i + 1] - _x[i];
                    VarType v = _f[i + 1] - _f[i];
                    _c2[i] = (3 * v - h * (2 * _c1[i] + _c1[i + 1]))
                             / (h * h);
                    _c3[i] = (-2 * v + h * (_c1[i] + _c1[i + 1]))
                             / (h * h * h);
                }
            }
        }
    }

    // Function evaluation code
    private
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
        this(VarArray x, FuncArray y, bool calcNow = true)
        {
            this(x.length);
            setAll(x, y);
        }
    }
} // TODO: unittest

/* -------------------------------------------------------------------------- */

// TODO: 1d B-spline
// TODO: non-spline interpolators
// TODO: 2d bilinear interpolation
// TODO: 2d bicubic interpolation

/* ========================================================================== */
/* NOTE: Perhaps these functions should go to a separte module */
private:

/* Common binary search in a sorted array.
   Returns index of the first element in the first interval that contains x */
private size_t binarySearch(A, X)(A a, X x)
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

private template ProduceArray(ElementOrContainer)
{
    static if(is(BaseElementType!ElementOrContainer == ElementOrContainer))
        alias ElementOrContainer[] ProduceArray;
    else
        alias ElementOrContainer ProduceArray;
}
