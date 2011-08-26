module scid.splines.univariate.base;

public import std.range;

public import scid.splines.common;
public import scid.splines.support;

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
