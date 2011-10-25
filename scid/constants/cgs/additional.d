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

// TODO: should this go to units modules?
/// Atomic mass constant [g]
enum double ATOMIC_MASS = 1.660_538_921e-24;
enum double D_ATOMIC_MASS = 0.000_000_073e-24; /// Standard uncertainty

/// Bohr magneton [erg / G]
enum double BOHR_MAGNETON = 9.274_009_68e-25;
enum double D_BOHR_MAGNETON = 0.000_000_20e-25; /// Standard uncertainty
