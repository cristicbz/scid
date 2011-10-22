/** Contains characteristics of elementary particles
  *
  * $(DDOC_SECTION_H System of units:) symmetric CGS
  * Origin:  $(LINK2 http://physics.nist.gov/cuu/index.html,
  *            The NIST Reference on Constants, Units, and Uncertainty)
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  * Version: alpha
  *
  */

module scid.constants.cgs.particles;

/// Electron characteristics
enum ELECTRON : double
{
    /// g-factor
    GFACTOR = -2.002_319_304_361_53,
    UC_GFACTOR = 0.000_000_000_000_53, /// Standard uncertainty

    /// Mass [g]
    MASS = 9.109_382_91e-28,
    UC_MASS = 0.000_000_40e-28 /// Standard uncertainty
}
