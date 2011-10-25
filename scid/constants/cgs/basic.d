/** Contains fundamental and other widely used constants
  *
  * $(DDOC_SECTION_H System of units:) symmetric CGS
  * Origin:  $(LINK2 http://physics.nist.gov/cuu/index.html,
  *            The NIST Reference on Constants, Units, and Uncertainty)
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  * Version: alpha
  *
  */

module scid.constants.cgs.basic;

/// Newtonian constant of gravitation, [cm^3 / (g * s^2)]
enum double GRAVITY = 6.673_84e-8;
enum double U_GRAVITY = 0.000_80e-8; /// Standard uncertainty

/// Speed of light in vacuum [cm / s]
enum double SPEEDLIGHT = 2.997_924_58e+10;

/// Planck constant [erg * s]
enum double PLANCK = 6.626_069_57e-27;
enum double D_PLANCK = 0.000_000_29e-27; /// Standard uncertainty

/// Planck constant over 2*pi [erg * s]
enum double PLANCK_2PI = 1.054_571_726e-27;
enum double D_PLANCK_2PI = 0.000_000_047e-27; /// Standard uncertainty

/// Elementary charge [statC]
enum double ECHARGE = 4.803_204_25e-10;
enum double D_ECHARGE = 0.000_000_10e-10; /// Standard uncertainty

/// Fine-structure constant
enum double FINE_STRUCT = 7.297_352_5698e-3;
enum double D_FINE_STRUCT = 0.000_000_0024e-3; /// Standard uncertainty

/// Avogadro constant
enum double AVOGADRO = 6.022_141_29e+23;
enum double D_AVOGADRO = 0.000_000_27e+23; /// Standard uncertainty

/// Boltzman constant [erg]
enum double BOLTZMAN = 1.380_6488e-16;
enum double D_BOLTZMAN = 0.000_0013e-16; /// Standard uncertainty
