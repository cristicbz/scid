/** Contains characteristics of elementary particles
  *
  * $(DDOC_SECTION_H System of units:) symmetric SI
  * Origin:  $(LINK2 http://physics.nist.gov/cuu/index.html,
  *            The NIST Reference on Constants, Units, and Uncertainty)
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  * Version: alpha
  *
  */

module scid.constants.si.particles;

/// Electron characteristics
enum ELECTRON : double
{
    /// Charge to mass quotient [C / kg]
    CHARGE_TO_MASS = -1.758_820_088e+11,
    UC_CHARGE_TO_MASS = 0.000_000_039e+11, /// Standard uncertainty

    /// Compton wavelength [m]
    COMPTONWL = 2.426_310_2389e-12,
    UC_COMPTONWL = 0.000_000_0016e-12, /// Standard uncertainty

    /// g-factor
    GFACTOR = -2.002_319_304_361_53,
    UC_GFACTOR = 0.000_000_000_000_53, /// Standard uncertainty

    /// Magnetic momentum [J / T]
    MAG_MOM = -9.284_764_30e-24,
    UC_MAG_MOM = 0.000_000_21e-24, /// Standard uncertainty

    /// Mass [kg]
    MASS = 9.109_382_91e-31,
    UC_MASS = 0.000_000_40e-31 /// Standard uncertainty
}

/// Proton characteristics
enum PROTON : double
{
    /// Charge to mass quotient [C / kg]
    CHARGE_TO_MASS = 9.578_833_58e+7,
    UC_CHARGE_TO_MASS = 0.000_000_21e+7, /// Standard uncertainty

    /// Compton wavelength [m]
    COMPTONWL = 1.321_409_856_23e-15,
    UC_COMPTONWL = 0.000_000_000_94e-15, /// Standard uncertainty

    /// g-factor
    GFACTOR = 5.585_694_713,
    UC_GFACTOR = 0.000_000_046, /// Standard uncertainty

    /// Magnetic momentum [J / T]
    MAG_MOM = 1.410_606_743e-26,
    UC_MAG_MOM = 0.000_000_033e-26, /// Standard uncertainty

    /// Mass [kg]
    MASS = 1.672_621_777e-27,
    UC_MASS = 0.000_000_074e-27 /// Standard uncertainty
}

/// Neutron characteristics
enum NEUTRON : double
{
    /// Compton wavelength [m]
    COMPTONWL = 1.319_590_9068e-15,
    UC_COMPTONWL = 0.000_000_0011e-15, /// Standard uncertainty

    /// g-factor
    GFACTOR = -3.826_085_45,
    UC_GFACTOR = 0.000_000_90, /// Standard uncertainty

    /// Magnetic momentum [J / T]
    MAG_MOM = -0.966_236_47e-26,
    UC_MAG_MOM = 0.000_000_23e-26, /// Standard uncertainty

    /// Mass [kg]
    MASS = 1.674_927_351e-27,
    UC_MASS = 0.000_000_074e-27 /// Standard uncertainty
}
