/** Contains fundamental world constants
  *
  * $(DDOC_SECTION System of units: symmetric CGS)
  * Origin:  $(LINK2 http://physics.nist.gov/cuu/index.html,
  *            The NIST Reference on Constants, Units, and Uncertainty)
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  * Version: alpha
  *
  */

module scid.constants.cgs.fundamental;

/// Newtonian constant of gravitation, [cm^3 / (g * s^2)]
enum double GRAVITY = 6.673_84e-8;
enum double U_GRAVITY = 0.000_80e-8; /// Standard uncertanity

/// Speed of light in vacuum [cm / s]
enum double SPEEDLIGHT = 2.997_924_58e+10;

enum PLANCK : double
{
    /// Planck constant [erg * s]
    H = 6.626_069_57e-27,
    U_H = 0.000_000_29e-27, /// Standard uncertanity

    /// Planck constant over 2*pi [erg * s]
    HBAR = 1.054_571_726e-27,
    U_HBAR = 0.000_000_047e-27 /// Standard uncertanity
}

/// Elementary charge [statC]
enum double ECHARGE = 4.803_204_25e-10;
enum double U_ECHARGE = 0.000_000_10e-10; /// Standard uncertanity
