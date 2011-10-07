/** Contains common features for 3rd-order splines.
  *
  * Version: 0.7-b
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev.
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  */
module scid.splines.univariate.poly3;

package mixin template poly3Base(VarType, FuncType)
{
    // Interpolant storage
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

        void _allocContents(size_t maxSize)
        {
            _c1.length = maxSize;
            _c2.length = maxSize - 1;
            _c3.length = maxSize - 1;
            static if(optim == SplineOptim.fixVar)
            {
                _allocOptFixVar(maxSize);
            }
            _maxSize = maxSize;
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

    public @property
    {
        /** Spline coefficients (read only)
          */
        const(FuncType)[] c1()
        {
            return _c1[0..(pointsNum - 1)];
        }

        ///ditto
        const(FuncType)[] c2()
        {
            return _c2[0..(pointsNum - 1)];
        }

        ///ditto
        const(FuncType)[] c3()
        {
            return _c3[0..(pointsNum - 1)];
        }
    }
}
