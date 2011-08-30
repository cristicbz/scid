/** Contains common features used by all splines.
  *
  * Version: 0.4-a
  * Authors: Maksim Zholudev
  * Copyright: Copyright (c) 2011, Maksim Zholudev.
  * License: $(WEB boost.org/LICENSE_1_0.txt, Boost License 1.0)
  */
module scid.splines.common;

/// The type of optimization
enum SplineOptim
{
    normal, /// No special features.
    fixVar /** Preliminary calculations depending only on VVA are made if it
             * changes. Minimal amount of operations is performed when FVA
             * changes.
             */
}
