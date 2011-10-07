/** Provides univariate linear spline.
  *
  * Version: 0.7-b
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev.
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  */
module scid.splines.univariate.linear;

import scid.splines.univariate.base;

/** Linear one-dimensional spline (order = 1, defect = 1).
  *
  * No boundary conditions are supported.
  *
  * Params:
  *     EocVar = type of variable or variable value array (VVA)
  *     EocFunc = type of function or function value array (FVA)
  *     optim = spline optimization type (normal by default)
  */
struct SplineLinear(EocVar, EocFunc,
                    SplineOptim optim = SplineOptim.normal)
{
    mixin splineBase!(EocVar, EocFunc);

    // Interpolant and optimization variables storage
    private
    {
        // Spline parameters:
        FuncType[] _c1;
        /* The interpolant is:
         *     f(x) = _f[i] + _c1[i] * dx
         *     dx = x - _x[i]
         */

        static if(optim == SplineOptim.fixVar)
        {
            /* Data depending only on variable values */
            private VarType[] _rdx;
        }

        void _allocContents(size_t maxSize)
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
    package
    {
        static if(optim == SplineOptim.fixVar)
        {
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
    package
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

        /** Spline coefficients (read only)
          */
        @property const(FuncType)[] c1()
        {
            return _c1[0..(pointsNum - 1)];
        }
    }
}

unittest
{
    double[] x = [1, 2, 3];
    double[] y = [-1, 1, -1];
    auto spl = SplineLinear!(double, double)(x, y);
    assert(spl._calcFunction(1.5, 0) == 0);
    assert(spl._calcDeriv(1.5, 0) == 2);
    assert(spl._calcFunction(2.75, 1) == -0.5);
    assert(spl._calcDeriv(2.75, 1) == -2);
}

unittest
{
    double[] x = [1, 2, 3];
    double[] y = [-1, 1, -1];
    auto spl = SplineLinear!(double, double, SplineOptim.fixVar)(x, y);
    assert(spl._calcFunction(1.5, 0) == 0);
    assert(spl._calcDeriv(1.5, 0) == 2);
    assert(spl._calcFunction(2.75, 1) == -0.5);
    assert(spl._calcDeriv(2.75, 1) == -2);
}
