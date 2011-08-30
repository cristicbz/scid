/** Provides univariate akima spline.
  *
  * Version: 0.4-a
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev.
  * License: $(WEB boost.org/LICENSE_1_0.txt, Boost License 1.0)
  */
module scid.splines.univariate.akima;

import std.math;

import scid.splines.univariate.base;

/** One-dimensional Akima interpolation.
  *
  * Params:
  *     EocVar = type of variable
  *     EocFunc = type of function
  *     optim = spline optimization type (normal by default)
  *     Props = strings, describing properties of FuncType type to process.
  *             It is necessary because for this kind of splines not only
  *             linear operations are performed with function values.
  */
  /* FIXME:
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
    // FIXME: Implement support of compound types
    static assert(Props.length == 0,
                  "Compound types support is not implemented yet");

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
                                  "Properties must be strings");
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
                /* TODO: After testing, refuse the workspaces by using
                 *       the coefficient arrays instead.
                 */
                size_t N = pointsNum - 1;
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
}

// TODO: unittest
