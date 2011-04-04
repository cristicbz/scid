/** This module contains classes for interpolation.
  *
  * This is a pilot version. Everything (even the interface) can be changed.
  *
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev. All rights reserved.
  * License: Boost License 1.0
  */

module scid.spline;

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

// After gsl_interp_accel
/** This structure helps to avoid repeating search.
  * It caches the previous value of an index lookup. When the subsequent
  * interpolation point falls in the same interval its index value can be
  * returned immediately.
  * In Spline class each method that evaluates function or derivative has two
  * versions: with accelerator and without it.
  */
/* TODO: This structure is supposed to be opaque for Spline classes.
 *       Should it go to a separate module?
 */
/* TODO: May be a better way is to provide a class which has an
 *       accelerator and a link to Spline object inside?
 */
struct SplineSearchAccel(Tvar)
{
  private size_t index = 0; // The index of current interval
  /** Find the value x in a sorted array a
    */
  size_t find(Tvar[] a, Tvar x) // TODO: Should it be private too?
  {
    // Does current interval contain x?
    if((a[index + 1] < x) || (a[index] > x))
      index = binarySearch(a, x); // Find new interval if it doesn't
    return index;
  }
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
  // TODO: add functions for derivatives and integrals
  abstract Tfunc _calcFunction(Tvar x, size_t index);
  // Minimal number of points needed to build the spline
  // TODO: Is it possible to make it variable not function?
  abstract size_t minPoints();
public:
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
    assert(f.length == x.length);
    assert(x.length >= minPoints);
    for(size_t i = 0; i < x.length-1; ++i)
      assert(x[i] < x[i+1]);
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
    assert(x >= _x[0]);
    assert(x <= _x[$-1]);
  }
  body
  {
    return _calcFunction(x, binarySearch(_x, x));
  }
  Tfunc eval(Tvar x, ref SplineSearchAccel!(Tvar) accel)
  in
  {
    assert(x >= _x[0]);
    assert(x <= _x[$-1]);
  }
  body
  {
    return _calcFunction(x, accel.find(_x, x));
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

  size_t minPoints()
  {
    return 2;
  }

public:
  this(Tfunc[] f, Tvar[] x, bool duplicate = false)
  {
    super(f, x, duplicate);
  }
}
