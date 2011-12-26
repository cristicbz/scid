module scid.demos.common;

version( demo ):
import std.math;
import std.string;

public import std.conv;
public import std.exception;
public import std.stdio;
public import std.typetuple;

public import scid.matvec;
public import scid.common.traits;
public import scid.common.meta;

/** myApproxEqual that works on complex elements as well */
bool myApproxEqual( T, U )( T a, U b ) {
    static if( is( T : cdouble ) && is( U : cdouble ) )
        return myApproxEqual( a.re, b.re ) && myApproxEqual( a.im, b.im );
    else
        return approxEqual( a, b );
}

/** Enforce the data & dimensions of a matrix. */
void enforceMatData( M, E )( auto ref M m, size_t r, size_t c, E[] expected ) {
    //debug {
        enum epsilon = 1e-3;
        enforce( m.rows == r, format("Wrong no. of rows %d vs %d", m.rows, r) );
        enforce( m.columns == c, format("Wrong no. of rows %d vs %d", m.columns, c ) );
        auto a = m.cdata[ 0 .. expected.length ];
        auto b = to!(BaseElementType!M[])(expected);
        foreach( i ; 0 .. a.length ) {
            enforce( myApproxEqual( a[i], b[i] ),
               "Expected " ~ to!string(expected) ~ ", got "~ to!string(a) ~
                " (" ~ to!string(i) ~  ")"
            );
        }
    //}
}

/** Enforce the data and length of a vector. */
void enforceVecData( V, E )( auto ref V v, E[] expected ) {
    //debug {
        enforce( v.length == expected.length, format("Wrong vector length: %d vs. %d", v.length, expected.length) );
        enforce( v.cdata[ 0 .. expected.length ] == to!(BaseElementType!V)(expected) );
    //}
}

/** Enforce that a matrix is myApproxEqual to a built-in matrix. */
void enforceMatApproxEqual( int line = __LINE__, A, E )( auto ref A a, E[][] b ) {
    string message() { return "At " ~ to!string(line) ~ ": Expected " ~ to!string(b) ~ ", got " ~ a.toString() ~ ". "; }

    enforce( a.rows == b.length, message() ~ "First dimensions do not match." );

    if( a.rows == 0 )
        return;

    enforce( a.columns == b[0].length, message() ~ "Second dimensions do not match." );

    foreach( i ; 0 .. a.rows ) {
        foreach( j ; 0 .. a.columns ) {
            enforce( myApproxEqual( a[ i, j ], b[ i ][ j ] ), message() ~
                format("Elements at position (%d,%d) do not match.", i, j)
            );
        }
    }
}

// Separating MatrixTypes into row and column major and doing the
// syntax test in separate files.  Otherwise the compiler runs
// out of memory on this test.
// Some of these types have to be disabled otherwise the compiler runs out of memory.
template MatrixTypesColumn( T ) {
    alias TypeTuple!(
         Matrix!T,
         TriangularMatrix!T,
         TriangularMatrix!(T, MatrixTriangle.Lower ),
         SymmetricMatrix!T,
         SymmetricMatrix!(T, MatrixTriangle.Lower ),
        //DiagonalMatrix!T // Issue 81
    ) MatrixTypesColumn;
}

template MatrixTypesRow( T ) {
    alias TypeTuple!(
         Matrix!(T,StorageOrder.RowMajor),
         TriangularMatrix!(T, MatrixTriangle.Upper, StorageOrder.RowMajor),
         TriangularMatrix!(T, MatrixTriangle.Lower, StorageOrder.RowMajor),
         SymmetricMatrix!(T, MatrixTriangle.Upper, StorageOrder.RowMajor),
         SymmetricMatrix!(T, MatrixTriangle.Lower, StorageOrder.RowMajor),
    ) MatrixTypesRow;
}

// Ideally we should test more types, but if we do the
// compiler runs out of memory.
alias TypeTuple!(double) ElementTypes;

/** Check that empty matrices and vectors don't cause crashes. */
void emptyMatVecTest()() {
    Matrix!double m1,m2;
    Vector!double v1,v2;

    // We're really checking there's no assertions or memory problems
    // Empty vectors
    v1[] = v2;
    v2[] = v1[ 0 .. 0 ];
    v1[] *= 2.0;
    v1[] += v2;
    double x = eval( v1.t*v2 );

    // Empty matrices
    m1[ 0 .. 0 ][ 0 .. 0 ] = m2;
    m1[] += m2;
    m1[] *= m2;

    // The inverse of an empty matrix is the empty matrix. After all, it satisfies
    // the properties of the identity matrix.
    enforce( eval( inv( m1 ) ).empty );
    enforce( ( eval( inv( m1 ) * m2 ) ).empty );
}
