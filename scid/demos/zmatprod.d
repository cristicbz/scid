module scid.demos.zmatprod;

version( demo ):
import scid.demos.common;

void zMatProdTest() {
    alias Matrix!cdouble                            ColMat;
    alias Matrix!( cdouble, StorageOrder.RowMajor ) RowMat;

    auto aInit    = [ [  1. + 1.i,  2. - 1.i ], [ -1. + 2.i,  1. + 1.i ] ];
    auto bInit    = [ [ -1. - 1.i,  2. - 2.i ], [  2. + 1.i,  0. - 1.i ] ];
    auto correct1 = [ [  5. - 2.i,  3. - 2.i ], [  4. + 2.i,  3. + 5.i ] ];
    auto correct2 = [ [  2. - 4.i,  3. - 4.i ], [  1. + 1.i,  6. + 1.i ] ];
    auto correct3 = [ [ -2. - 5.i, -2. - 3.i ], [  2. - 4.i,  5. - 3.i ] ];
    auto correct4 = [ [  4. + 2.i,  4. + 3.i ], [ -1. + 1.i, -1. + 6.i ] ];

    auto colA = ColMat( aInit ), colB = ColMat( bInit );
    auto rowA = RowMat( aInit ), rowB = RowMat( bInit );
    ColMat colC;
    RowMat rowC;

    // check all cases
    colC[] = colA * colB;     enforceMatApproxEqual( colC, correct1 );
    colC[] = colA * rowB;     enforceMatApproxEqual( colC, correct1 );
    colC[] = rowA * colB;     enforceMatApproxEqual( colC, correct1 );
    colC[] = rowA * rowB;     enforceMatApproxEqual( colC, correct1 );

    colC[] = colA.t * colB.t; enforceMatApproxEqual( colC, correct2 );
    colC[] = colA.t * rowB.t; enforceMatApproxEqual( colC, correct2 );
    colC[] = rowA.t * colB.t; enforceMatApproxEqual( colC, correct2 );
    colC[] = rowA.t * rowB.t; enforceMatApproxEqual( colC, correct2 );

    rowC[] = colA * colB;     enforceMatApproxEqual( rowC, correct1 );
    rowC[] = colA * rowB;     enforceMatApproxEqual( rowC, correct1 );
    rowC[] = rowA * colB;     enforceMatApproxEqual( rowC, correct1 );
    rowC[] = rowA * rowB;     enforceMatApproxEqual( rowC, correct1 );

    rowC[] = colA.t * colB.t; enforceMatApproxEqual( rowC, correct2 );
    rowC[] = colA.t * rowB.t; enforceMatApproxEqual( rowC, correct2 );
    rowC[] = rowA.t * colB.t; enforceMatApproxEqual( rowC, correct2 );
    rowC[] = rowA.t * rowB.t; enforceMatApproxEqual( rowC, correct2 );

    colC[] = colA.t * colB;   enforceMatApproxEqual( colC, correct3 );
    colC[] = colA.t * rowB;   enforceMatApproxEqual( colC, correct3 );
    colC[] = rowA.t * colB;   enforceMatApproxEqual( colC, correct3 );
    colC[] = rowA.t * rowB;   enforceMatApproxEqual( colC, correct3 );

    colC[] = colA * colB.t;   enforceMatApproxEqual( colC, correct4 );
    colC[] = colA * rowB.t;   enforceMatApproxEqual( colC, correct4 );
    colC[] = rowA * colB.t;   enforceMatApproxEqual( colC, correct4 );
    colC[] = rowA * rowB.t;   enforceMatApproxEqual( colC, correct4 );

    rowC[] = colA.t * colB;   enforceMatApproxEqual( rowC, correct3 );
    rowC[] = colA.t * rowB;   enforceMatApproxEqual( rowC, correct3 );
    rowC[] = rowA.t * colB;   enforceMatApproxEqual( rowC, correct3 );
    rowC[] = rowA.t * rowB;   enforceMatApproxEqual( rowC, correct3 );

    rowC[] = colA * colB.t;   enforceMatApproxEqual( rowC, correct4 );
    rowC[] = colA * rowB.t;   enforceMatApproxEqual( rowC, correct4 );
    rowC[] = rowA * colB.t;   enforceMatApproxEqual( rowC, correct4 );
    rowC[] = rowA * rowB.t;   enforceMatApproxEqual( rowC, correct4 );
}
