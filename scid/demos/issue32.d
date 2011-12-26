module scid.demos.issue32;

version( demo ):
import scid.demos.common;

/** Issue 32 - D'tor problems/assert failures, Linux Only */
void testIssue32()() {
    auto xTx = Matrix!double(1, 1);
    auto xTy = Matrix!double(1, 2);
    xTx[0, 0] = 31;

    xTy[0, 0] = 41;
    xTy[0, 1] = 59;

    auto ret = eval(inv(xTx) * xTy);
}
