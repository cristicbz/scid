/** This module contains types for interpolation.
  *
  * This is an alpha version. The interface can be changed.
  *
  * Version: 0.1a
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev. All rights reserved.
  * License: Boost License 1.0
  */
module scid.interpolation;
// TODO: add more embedded documentation
// FIXME: add "in" and "out" qualifiers for function arguments

// Common binary search in a sorted array
// FIXME: Does phobos already have such function?
private pure size_t BinarySearch(T)(T[] a, T x)
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

/** Structure to wrap a one-dimensional spline.
  * It caches the previous value of an index lookup. When the subsequent
  * interpolation point falls in the same interval its index value can be
  * returned immediately.
  * Many Spline1dView structures can refer the same spline.
  * Params:
  *   Tspline = type of spline to wrap
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
      _index = BinarySearch(_spline._x, x); // Find new interval if it doesn't
  }

public:
  this(ref Tspline spline)
  {
    _spline = &spline;
  }

  /** Evaluate function value at a given point
    * Params:
    *   x = variable value
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

/** Linear one-dimensional spline (order = 1, defect = 1).
  * Params:
  *   Tfunc = type of function
  *    Tvar = type of variable
  */
struct SplineLinear(Tfunc, Tvar)
{
private:
  Tfunc[] _f; // function value array
  Tvar[] _x; // variable value array

  // Spline parameters:
  Tfunc[] _c1;
  /* The interpolant is:
   *   f(x) = _f[i] + _c1[i] * dx
   *   dx = x - _x[i]
   */

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
  /// Flags used by spline creation function
  enum flags
  {
    dupSource = 0x00_00_00_01, // duplicate point arrays
    update = 0x00_00_00_02 // calculate spline coefficients immediately
  }
  /// Minimal number of points needed for the spline
  enum size_t minPoints = 2;

  this(Tfunc[] f, Tvar[] x, flags flg = flags.update)
  {
    bind(f, x, flg);
  }

  /** Select point arrays
    * Params:
    *   f = array of function values
    *   x = array of variable values
    *       must be sorted in ascending order if update flag is on
    */
  void bind(Tfunc[] f, Tvar[] x, flags flg = flags.update)
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
      update();
  }

  /** Calculate spline prarmeters
    */
  void update()
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
  body
  {
    _c1.length = _x.length-1;
    for(size_t i = 0; i < _x.length-1; ++i)
      _c1[i] = (_f[i+1] - _f[i]) / (_x[i+1] - _x[i]);
  }
}

// TODO: use SciD functions instead
// Solve Ax = b linear system for three-diagonal matrix A
// NOTE: x is passed by refence to allow this function allocate this array
private void Sweep3Diag(T)(in T[][3] a, in T[] b, T[] x)
in
{
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

// TODO: cyclic sweep for periodic splines

/** Cubic one-dimensional spline (order = 3, defect = 1).
  * Params:
  *   Tfunc = type of function
  *    Tvar = type of variable
  */
// NOTE: not tested properly
// FIXME: boundary conditions
struct SplineCubic(Tfunc, Tvar)
{
private:
  Tfunc[] _f; // function value array
  Tvar[] _x; // variable value array

  // Spline parameters:
  Tfunc[] _c1;
  Tfunc[] _c2;
  Tfunc[] _c3;
  /* The interpolant is:
   *   f(x) = _f[i] + _c1[i] * dx + _c2[i] * dx*dx + _c3[i] * dx*dx*dx
   *   dx = x - _x[i]
   */

  // Calculate function in a given interval
  Tfunc _calcFunction(Tvar x, size_t index)
  {
    double dx = x - _x[index];
    return _f[index] + dx * (_c1[index] + dx * (_c2[index] + dx * _c3[index]));
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
    // ignore boundary conditions and build a periodic spline
    periodic = 0x00_00_00_03
  }
  /// Boundary conditions by spline creation function
  // FIXME: not implemented so far
  enum boundCond
  {
    deriv1,
    deriv2
  }
  /// Minimal number of points needed for the spline
  enum size_t minPoints = 3;

  this(Tfunc[] f, Tvar[] x, flags flg = flags.update)
  {
    bind(f, x, flg);
  }

  /** Select point arrays
    * Params:
    *   f = array of function values
    *   x = array of variable values
    *       must be sorted in ascending order if update flag is on
    */
  void bind(Tfunc[] f, Tvar[] x, flags flg = flags.update)
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
      update();
  }

  /** Calculate spline prarmeters
    */
  void update()
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

    // FIXME: use boundary conditions here
    a[1][0] = 1;
    a[2][0] = 0;
    b[0] = 0;

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

    // FIXME: use boundary conditions here
    a[0][$-1] = 0;
    a[1][$-1] = 1;
    b[$-1] = 0;

    _c1.length = _x.length;
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

// TODO: 1d Akima spline
// TODO: 1d B-spline
// TODO: non-spline interpolators
// TODO: 2d bilinear interpolation
// TODO: 2d bicubic interpolation
