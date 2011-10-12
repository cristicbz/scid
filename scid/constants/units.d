/** Contains conversion coefficients for units.
  *
  * Origin: * $(LINK2 http://physics.nist.gov/cuu/index.html,
                The NIST Reference on Constants, Units, and Uncertainty)
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  *
  */

module scid.constants.units;

///???

//FIXME
    /*
     *
     */
enum double nano = 1e-9;
    /*
     *
     */

//------------------------------------------------------------------------------

/// Length units
enum cm : double
{
    Angstrom = 1e8
}

/// ditto
enum Angstrom : double
{
    cm = 1e-8, ///centimeter
    m = 1e-10, ///meter
}

//------------------------------------------------------------------------------

/// Energy units
enum erg : double
{
    eV = 1.0, // FIXME
    meV = 1.0 // FIXME
}
