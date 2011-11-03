/** Contains SI prefixes for multiple and submultiple units.
  *
  * Origin: Wikipedia
  *
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  * Version: alpha
  *
  * Exmaples:
  * ----------
  * TimeInPicoseconds = TimeInSeconds / PREFIX.PICO;
  * ----------
  *
  */

module scid.constants.units.siprefixes;

/// SI prefixes
enum PREFIX : double
{
    YOTTA = 1e24, /// Multiple prefixes
    ZETTA = 1e21, /// ditto
    EXA   = 1e18, /// ditto
    PETA  = 1e15, /// ditto
    TERA  = 1e12, /// ditto
    GIGA  = 1e9, /// ditto
    MEGA  = 1e6, /// ditto
    KILO  = 1e3, /// ditto
    HECTO = 1e2, /// ditto
    DECA  = 1e1, /// ditto

    DECI  = 1e-1, /// Submultiple prefixes
    CENTI = 1e-2, /// ditto
    MILLI = 1e-3, /// ditto
    MICRO = 1e-6, /// ditto
    NANO  = 1e-9, /// ditto
    PICO  = 1e-12, /// ditto
    FEMTO = 1e-15, /// ditto
    ATTO  = 1e-18, /// ditto
    ZEPTO = 1e-21, /// ditto
    YOCTO = 1e-24 /// ditto
}
