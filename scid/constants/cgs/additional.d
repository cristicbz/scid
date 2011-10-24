/** Contains some useful constants
  *
  * $(DDOC_SECTION_H System of units:) symmetric CGS
  * Origin:  $(LINK2 http://physics.nist.gov/cuu/index.html,
  *            The NIST Reference on Constants, Units, and Uncertainty)
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  * Version: alpha
  *
  */

module scid.constants.cgs.additional;

/// Atomic mass constant [g]
enum double ATOMIC_MASS = 1.660_538_921e-24;
enum double D_ATOMIC_MASS = 0.000_000_073e-24; /// Standard uncertainty

/// Avogadro constant
enum double AVOGADRO = 6.022_141_29e+23;
enum double D_AVOGADRO = 0.000_000_27e+23; /// Standard uncertainty

/// Bohr magneton [erg / G]
enum double BOHR_MAGNETON = 927.400_968e-23;
enum double D_BOHR_MAGNETON = 0.000_020e-23; /// Standard uncertainty

/// Boltzman constant [erg]
enum double BOLTZMAN = 1.380_6488e-16;
enum double D_BOLTZMAN = 0.000_0013e-16; /// Standard uncertainty
