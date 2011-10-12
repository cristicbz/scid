/** Contains fundamental world constants in symmetric CGS system
  *
  * Origin: * $(LINK2 http://physics.nist.gov/cuu/index.html,
                The NIST Reference on Constants, Units, and Uncertainty)
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  *
  */

/// Newtonian constant of gravitation, [cm^3 / (g * s^2)]
enum double GravConst = 6.673_84e-8;
enum double u_GravConst = 0.000_80e-8; /// Standard uncertanity

/// Speed of light in vacuum [cm / s]
enum double SpeedLight = 2.997_924_58e+10;

enum PlanckConst : double
{
    /// Planck constant [erg * s]
    h = 6.626_069_57e-27;
    u_h = 0.000_000_29e-27; /// Standard uncertanity

    /// Planck constant over 2*pi [erg * s]
    hbar = 1.054_571_726e-27;
    u_hbar = 0.000_000_047e-27; /// Standard uncertanity
}

/// Elementary charge [statC]
enum double ECharge = 4.803_204_25e-10;
enum double u_ECharge = 0.000_000_10e-10; /// Standard uncertanity
