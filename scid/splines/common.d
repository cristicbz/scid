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
