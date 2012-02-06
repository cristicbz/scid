module scid.demos.zmatprod;

version( demo ):
import scid.demos.common;

void zMatProdTest() {
    alias Matrix!cdouble                            ColMat;
    alias Matrix!( cdouble, StorageOrder.RowMajor ) RowMat;

    auto aInit    = [ [  1. + 1i,  2. - 1i ], [ -1. + 2i,  1. + 1i ] ];
    auto bInit    = [ [ -1. - 1i,  2. - 2i ], [  2. + 1i,  0. - 1i ] ];
    auto correct1 = [ [  5. - 2i,  3. - 2i ], [  4. + 2i,  3. + 5i ] ];
    auto correct2 = [ [  2. - 4i,  3. - 4i ], [  1. + 1i,  6. + 1i ] ];
    auto correct3 = [ [ -2. - 5i, -2. - 3i ], [  2. - 4i,  5. - 3i ] ];
    auto correct4 = [ [  4. + 2i,  4. + 3i ], [ -1. + 1i, -1. + 6i ] ];

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
