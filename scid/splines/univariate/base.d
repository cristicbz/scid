module scid.splines.univariate.base;

public import std.range;

public import scid.splines.base;

/** Common functions for univariate splines.
  * "Eoc" prefix means Element Or Container
  */
package mixin template splineBase(EocVar, EocFunc)
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
    static assert(isRandomAccessRange!(VarArray), "Illegal variable type");
    static assert(isRandomAccessRange!(FuncArray), "Illegal function type");

    alias BaseElementType!VarArray VarType;
    alias BaseElementType!FuncArray FuncType;

    package
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
        /** Number of points.
          */
        @property const pointsNum()
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

        /** Test whether given point is inside spline range
          */
        bool pointInsideRange(VarType x)
        {
            /* TODO: change pointsNum to $
               when opDollar will be implemented */
            return (x >= _x[0]) && (x <= _x[pointsNum-1]);
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
                       "Variable value array must be strictly sorted ascending");
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
