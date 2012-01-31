/**
This module contains demos for SciD.  These also serve as integration
tests.  Most of the tests are templated so that they aren't compiled
unless called from main().  This is necessary to work around OptLink
bugs by skipping compiling/linking certain tests on Windows.
*/

version( demo ):

import std.stdio;
import scid.matvec;
import scid.storage.diagonalmat;
import scid.demos.basicexp;
import scid.demos.common;
import scid.demos.datainterface;
import scid.demos.dmatinv;
import scid.demos.dmatops;
import scid.demos.dmatprod;
import scid.demos.external;
import scid.demos.issue32;
import scid.demos.issue35;
import scid.demos.issue48;
import scid.demos.issue49;
import scid.demos.issue50;
import scid.demos.issue51;
import scid.demos.issue52;
import scid.demos.issue54;
import scid.demos.issue61;
import scid.demos.rangeinterface;
import scid.demos.zmat;
import scid.demos.zmatprod;



version(none) {
    // We can't systematically run all the tests on Windows because SciD is so
    // good at blowing up OptLink.
    void main() {
        auto x = DiagonalMatrix!double([1.,2,3,4]), y=DiagonalMatrix!double([2.,2.,2.,2.]);
        //  auto w = SymmetricMatrix!double([[1.,2.,3.,4],[1.,2.,3.,4],[1.,2.,3.,4],[1.,2.,3.,4]]);
        eval( x[0..2][] * y[][0..2] );
        //writeln(z.pretty);

        dMatInvTest();
        testIssue51();
        readln();
    }
} else {
    // Systematically run all of the demo/tests in this module.
    void main() {
         stderr.writeln("Starting demos...");
         basicExpressions();
         rangeInterface();
         dataInterface();
         externalViews();
         emptyMatVecTest();
         dMatOpsTest();
         dMatInvTest();
         dMatProdTest();
         zMatProdTest();
         zMatOpsTest();
         testIssue52();
         testIssue51();
         testIssue50();
         testIssue49();
         testIssue35();
         testIssue32();
         testIssue48();
//         testIssue54();  // This issue isn't fixed yet.
         testIssue61();
    }
}


