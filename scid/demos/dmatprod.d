module scid.demos.dmatprod;

version( demo ):
import scid.demos.common;

/** Specifically test all cases of matrix products for doubles (32 cases...). */
void dMatProdTest() {
    auto aInit    = [ [ 1.,  2, 1 ], [ 3., -1, 2]  ];
    auto bInit    = [ [ 1., -2 ],    [ -3., 2 ], [ -2., 1. ] ];
    auto correct1 = [ [-7., 3], [2., -6] ];
    auto correct2 = [ [-5.,3,1],[4.,-8,-5],[-3.,1,0] ];

    alias Matrix!double                            ColMat;
    alias Matrix!( double, StorageOrder.RowMajor ) RowMat;

    // check all cases except mixed transpose (dimensions wouldn't match)
    auto colA = ColMat( aInit ), colB = ColMat( bInit );
    auto rowA = RowMat( aInit ), rowB = RowMat( bInit );
    ColMat colC;
    RowMat rowC;

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

    // check mixed transpose
    auto correct3 = [[-8.,4],[5.,-6]];
    auto correct4 = [[-3.,1],[5.,-11]];
    colA[] = colA[0..2][0..2]; colB[] = colB[0..2][0..2];
    rowA[] = rowA[0..2][0..2]; rowB[] = rowB[0..2][0..2];

    colC[] = colA.t * colB; enforceMatApproxEqual( colC, correct3 );
    colC[] = colA.t * rowB; enforceMatApproxEqual( colC, correct3 );
    colC[] = rowA.t * colB; enforceMatApproxEqual( colC, correct3 );
    colC[] = rowA.t * rowB; enforceMatApproxEqual( colC, correct3 );

    colC[] = colA * colB.t; enforceMatApproxEqual( colC, correct4 );
    colC[] = colA * rowB.t; enforceMatApproxEqual( colC, correct4 );
    colC[] = rowA * colB.t; enforceMatApproxEqual( colC, correct4 );
    colC[] = rowA * rowB.t; enforceMatApproxEqual( colC, correct4 );

    rowC[] = colA.t * colB; enforceMatApproxEqual( rowC, correct3 );
    rowC[] = colA.t * rowB; enforceMatApproxEqual( rowC, correct3 );
    rowC[] = rowA.t * colB; enforceMatApproxEqual( rowC, correct3 );
    rowC[] = rowA.t * rowB; enforceMatApproxEqual( rowC, correct3 );

    rowC[] = colA * colB.t; enforceMatApproxEqual( rowC, correct4 );
    rowC[] = colA * rowB.t; enforceMatApproxEqual( rowC, correct4 );
    rowC[] = rowA * colB.t; enforceMatApproxEqual( rowC, correct4 );
    rowC[] = rowA * rowB.t; enforceMatApproxEqual( rowC, correct4 );
}

