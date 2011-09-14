/** Contains links to interpolation libraries.
  *
  * Version: 0.7-b
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev.
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
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
module scid.interpolation;

public import scid.splines.all;
