/** This module contains classes for interpolation.
  *
  * This is a pilot version. Everything (even the interface) can be changed.
  * This version is an alternative to the previous one in spline_a.d
  *
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev. All rights reserved.
  * License: Boost License 1.0
  */

module scid.spline_b;

// Common binary search in a sorted array
// TODO: Does phobos already have such function?
private size_t binarySearch(Tvar)(Tvar[] a, Tvar x)
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
  return ilo;
}

/** The base class for all one-dimensional splines
  * Tfunc - type of function
  * Tvar - type of variable
  */
class Spline(Tfunc, Tvar)
{
protected:
  Tfunc[] _f;
  Tvar[] _x;

  // Calculate the parameters for spline
  abstract void _calcParams();
  // Interpolation functions
  abstract Tfunc _calcFunction(Tvar x, size_t index);
  // TODO: add functions for derivatives and integrals
public:
  /** Minimal number of points needed to build the spline
    */
  static size_t minPoints;
  /** Calculate the parameters of the spline for given points
    * f[] - function values
    * x[] - points in ascending order
    * duplicate - whether to make a copy of arrays or not
    * If the arrays are not duplicated they should exist and not change
    * while the spline is used.
    */
  this(Tfunc[] f, Tvar[] x, bool duplicate = false)
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
    if(duplicate) // Make a copy of source data if required
    {
      _f = f.dup;
      _x = x.dup;
    }
    else // Otherwise just save links and hope it won't change
    {
      _f = f;
      _x = x;
    }
    _calcParams(); // Calculate spline paramwters
  }
  /** Evaluate the function value at the given point
    */
  Tfunc eval(Tvar x)
  in
  {
    assert((x >= _x[0]) && (x <= _x[$-1]),
           "Variable value out of interpolation range");
  }
  body
  {
    return _calcFunction(x, binarySearch(_x, x));
  }
}

// After gsl_interp_accel
/** This structure provides accelrated access to a spline.
  * It caches the previous value of an index lookup. When the subsequent
  * interpolation point falls in the same interval its index value can be
  * returned immediately.
  * Many accelerators can refer the same spline.
  */
struct SplineAccel(Tfunc, Tvar)
{
private:
  Spline!(Tfunc, Tvar) _spline; // Pointer to a spline object
  size_t _index = 0; // The index of current interval

  // Find the value x in a sorted array a
  void update(Tvar x)
  {
    // Does current interval contain x?
    if((_spline._x[_index + 1] < x) || (_spline._x[_index] > x))
      _index = binarySearch(_spline._x, x); // Find new interval if it doesn't
  }

public:
  this(Spline!(Tfunc, Tvar) spline)
  {
    _spline = spline;
  }

  Tfunc eval(Tvar x)
  in
  {
    assert((x >= _spline._x[0]) && (x <= _spline._x[$-1]),
           "Variable value out of interpolation range");
  }
  body
  {
    update(x);
    return _spline._calcFunction(x, _index);
  }
}

/** Linear interpolation
  */
class SplineLinear(Tfunc, Tvar) : Spline!(Tfunc, Tvar)
{
protected:
  Tfunc[] _c1; // the coefficients before x

  // Calculate the parameters for spline
  void _calcParams()
  {
    _c1.length = _x.length-1;
    for(size_t i = 0; i < _x.length-1; ++i)
      _c1[i] = (_f[i+1] - _f[i]) / (_x[i+1] - _x[i]);
  }

  // Interpolation functions
  Tfunc _calcFunction(Tvar x, size_t index)
  {
    return _f[index] + _c1[index] * (x - _x[index]);
  }
  // TODO: add functions for derivatives and integrals
public:
  static this()
  {
    minPoints = 2;
  }

  this(Tfunc[] f, Tvar[] x, bool duplicate = false)
  {
    super(f, x, duplicate);
  }
}
