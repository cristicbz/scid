/** Contains frequently used mathematical constants.
  * There are constants imported directly from std.math and some of them are
  * renamed.
  *
  * License: $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
  * Version: alpha
  *
  * Macros:
  *     SQRT = &radic;
  */

module scid.constants.math;

public import std.math :
    EXP1 = E,
    LOG2T,
    LOG2E,
    LOG2,
    LOG10E,
    LN2,
    LN10,
    PI,
    PI_2,
    PI_4,
    M_1_PI,
    M_2_PI,
    M_2_SQRTPI,
    SQRT2,
    SQRT1_2;

/* Values obtained from Wolfram Alpha.
 * Wolfram Alpha LLC. 2011. Wolfram|Alpha.
 * http://www.wolframalpha.com/input/?i=e+in+base+16
 * (access Nov 3, 2011).
 */
/** Square roots of small numbers
  */
enum SQRT : real
{
    _2 =  SQRT2,                                 /// $(SQRT)2
    _3 =  0x1.bb67ae8584caa73b25742d7078b84p+0L, /// $(SQRT)3
    _5 =  0x1.1e3779b97f4a7c15f39cc0605cedcp+1L, /// $(SQRT)5
    _6 =  0x1.3988e1409212e7d0321914321a556p+1L, /// $(SQRT)6
    _7 =  0x1.52a7fa9d2f8e9b78e753f30fe1bd1p+1L, /// $(SQRT)7
    _10 = 0x1.94c583ada5b529204a2bc830cd9c0p+1L, /// $(SQRT)10

    _1_2 =  SQRT1_2,                               /// $(SQRT)1/2
    _1_3 =  0x1.279a74590331c4d218f81e4afb258p-1L, /// $(SQRT)1/3
    _1_5 =  0x1.c9f25c5bfedd93565294670094afap-2L, /// $(SQRT)1/5
    _1_6 =  0x1.a20bd700c2c3dfc042cc1aed7871ep-2L, /// $(SQRT)1/6
    _1_7 =  0x1.83091e6a7f7e688a2cf23a5b4b213p-2L, /// $(SQRT)1/7
    _1_10 = 0x1.43d136248490edb36e896cf3d7b00p-2L  /// $(SQRT)1/10
}
