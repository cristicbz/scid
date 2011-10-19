/** Contains conversion coefficients for units.
  *
  * Origin:  $(LINK2 http://physics.nist.gov/cuu/index.html,
  *            The NIST Reference on Constants, Units, and Uncertainty)
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  * Version: alpha
  *
  */

module scid.constants.units;

/// SI prefixes

enum PREFIX : double
{
    YOTTA = 1e24,
    ZETTA = 1e21,
    EXA   = 1e18,
    PETA  = 1e15,
    TERA  = 1e12,
    GIGA  = 1e9,
    MEGA  = 1e6,
    KILO  = 1e3,
    HECTO = 1e2,
    DECA  = 1e1,

    DECI  = 1e-1,
    CENTI = 1e-2,
    MILLI = 1e-3,
    MICRO = 1e-6,
    NANO  = 1e-9,
    PICO  = 1e-12,
    FEMTO = 1e-15,
    ATTO  = 1e-18,
    ZEPTO = 1e-21,
    YOCTO = 1e-24
}

//------------------------------------------------------------------------------

/// Length units
enum CENTIMETER : double
{
    ANGSTROM = 1e8
}

/// ditto
enum ANGSTROM : double
{
    CENTIMETER = 1e-8, ///centimeter
    METER = 1e-10, ///meter
}

//------------------------------------------------------------------------------

/// Energy units
enum ERG : double
{
    ELECTRONVOLT = 1.0, // FIXME
    MILLIELECTRONVOLT = 1.0 // FIXME
}
