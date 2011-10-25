/** Contains conversion coefficients for units.
  *
  * Origin:  $(LINK2 http://physics.nist.gov/cuu/index.html,
  *            The NIST Reference on Constants, Units, and Uncertainty),
  *          Wikipedia
  *
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  * Version: alpha
  *
  * Exmaples:
  * ----------
  * LengthInMeters = LengthInAngstroms * ANGSTROM.METER;
  * ----------
  *
  * _PS postfix means "per second" ($BR)
  * _PM postfix means "per meter"
  */

module scid.constants.units.basic;

private
{
    import std.math;
    import scid.constants.cgs.basic;
}

//------------------------------------------------------------------------------
// Length units
public
{

/** Standard lenght unit in SI.
  *
  * Conversion coefficients are available for:
  */
enum METER : double
{
    ANGSTROM = 1e10 ///
}

/** Standard lenght unit in CGS.
  *
  * Conversion coefficients are available for:
  */
enum CENTIMETER : double
{
    ANGSTROM = 1e8 ///
}

/** Lenght unit.
  *
  * Conversion coefficients are available for:
  */
enum ANGSTROM : double
{
    METER = 1e-10, ///
    CENTIMETER = 1e-8 /// ditto
}

}

//------------------------------------------------------------------------------
// Energy units
public
{

/** Standard energy unit in SI.
  *
  * Conversion coefficients are available for:
  */
enum JOULE : double
{
    ERG = 1e7, ///
    ELECTRONVOLT = 6.241_509_343e+18, /// ditto
    MILLIELECTRONVOLT = 6.241_509_343e+15 /// ditto
}

/** Standard energy unit in CGS.
  *
  * Conversion coefficients are available for:
  */
enum ERG : double
{
    JOULE = 1e-7, ///
    ELECTRONVOLT = 6.241_509_343e+11, /// ditto
    MILLIELECTRONVOLT = 6.241_509_343e+8 /// ditto
}

/** Energy unit.
  *
  * Conversion coefficients are available for:
  */
enum ELECTRONVOLT : double
{
    JOULE = 1.602_176_565e-19, ///
    ERG = 1.602_176_565e-12, /// ditto
}

}

//------------------------------------------------------------------------------
// Electric charge units
public
{

/** Standard electric charge unit in SI.
  *
  * Conversion coefficients are available for:
  */
enum COULOMB : double
{
    STATCOULOMB = SPEEDLIGHT * 1e-1 ///
}

/** Standard electric charge unit in CGS.
  *
  * Conversion coefficients are available for:
  */
enum STATCOULOMB : double
{
    COULOMB = 1e1 / SPEEDLIGHT ///
}

}

//------------------------------------------------------------------------------
// Electric current units
public
{

/** Standard electric current unit in SI.
  *
  * Conversion coefficients are available for:
  */
enum AMPERE : double
{
    STATC_PS = COULOMB.STATCOULOMB ///
}

/** Standard electric current unit in CGS.
  *
  * Conversion coefficients are available for:
  */
enum STATC_PS : double
{
    AMPERE = STATCOULOMB.COULOMB ///
}

}

//------------------------------------------------------------------------------
// Electric voltage units
public
{

/** Standard electric voltage unit in SI.
  *
  * Conversion coefficients are available for:
  */
enum VOLT : double
{
    STATVOLT = STATCOULOMB.COULOMB * 1e9 ///
}

/** Standard electric voltage unit in CGS.
  *
  * Conversion coefficients are available for:
  */
enum STATVOLT : double
{
    VOLT = COULOMB.STATCOULOMB * 1e-9 ///
}

}

//------------------------------------------------------------------------------
// Magnetic induction units
public
{

/** Standard magnetic induction unit in SI.
  *
  * Conversion coefficients are available for:
  */
enum TESLA : double
{
    GAUSS = 1e4 ///
}

/** Standard magnetic induction unit in CGS.
  *
  * Conversion coefficients are available for:
  */
enum GAUSS : double
{
    TESLA = 1e-4 ///
}

}

//------------------------------------------------------------------------------
// Magnetic field strength units
public
{

/** Standard magnetic field strength unit in SI.
  *
  * Conversion coefficients are available for:
  */
enum AMPER_PM : double
{
    OERSTED = PI * 4e3 ///
}

/** Standard magnetic field strength unit in CGS.
  *
  * Conversion coefficients are available for:
  */
enum OERSTED : double
{
    AMPER_PM = 2.5e-4 / PI ///
}

}
