/** This module contains the MatrixView class as well as some
    functions related to it.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.matrix;


private import std.conv;
import std.array : appender;  // For toString
import std.exception;
import std.string: format, repeat;
import std.traits;
import std.typetuple;

import scid.core.meta;
import scid.core.traits;

import scid.matrixops;

version(unittest) {
    import scid.core.testing;
    import std.math;
}

/** Various matrix representations. */
enum Storage
{
    general,        /// General (dense) matrices
    triangular,     /// Packed storage of triangular matrices
    symmetric,      /// Packed storage of symmetric matrices
    diagonal,        /// Packed storage of diagonal matrices

    // These guys are kept for now for backwards compatibility but scheduled
    // for deprecation/removal.  The new convention is that enum members are
    // lowerCamelCase.

    General = general,
    Triangular = triangular,
    Symmetric = symmetric
}


/** In packed storage (triangular, symmetric, and Hermitian matrices),
    one can choose to store either the upper or lower triangle.
*/
enum Triangle : char
{
    upper = 'U',    /// Store upper triangle
    lower = 'L',     /// Store lower triangle

    // Scheduled for deprecation/removal
    Upper = upper,
    Lower = lower
}


/** A convenience function that allocates heap memory for a matrix,
    optionally sets the values of the matrix elements, and returns
    a MatrixView of the allocated memory.

    Examples:
    ---
    // Allocate general dense 3x4 matrix:
    auto denseMatrix = matrix!real(3, 4);

    // Allocate dense 3x2 zero-filled matrix:
    auto denseZeroMatrix = matrix!real(3, 2, 0.0L);

    // Allocate lower triangular 3x3 matrix:
    auto loMatrix = matrix!(real, Storage.triangular, Triangle.lower)(3);

    // Allocate upper triangular 2x2 matrix where the upper
    // triangular elements are set to 3.14.
    auto upMatrix = matrix!(real, Storage.triangular)(2, 3.14L);
    ---
*/
MatrixView!(T) matrix (T) (size_t rows, size_t cols) pure
{
    return typeof(return)(new T[rows*cols], rows, cols);
}


/// ditto
MatrixView!(T) matrix(T) (size_t rows, size_t cols, T init) pure
{
    auto array = new T[rows*cols];
    array[] = init;
    return typeof(return)(array, rows, cols);
}


///ditto
MatrixView!(T, stor, tri) matrix
    (T, Storage stor, Triangle tri = Triangle.upper)
    (size_t n, T init=T.init)
    pure
    if (stor == Storage.triangular)
{
    auto array = new T[(n*n+n)/2];
    if (init != T.init) array[] = init; // Because of DMD bug #3576 this can't
                                        // be done with a function overload.
    return typeof(return)(array, n);
}

/**
Convenience function for matrix literals or creating a matrix from an array
of arrays.
*/
MatrixView!(T) matrix(T)(const T[][] arrayOfArrays)
{
    typeof(return) ret;
    alias arrayOfArrays aa;  // Save typing.
    if(aa.length == 0) return ret;

    ret.array.length = aa.length * aa[0].length;
    ret.rows = aa.length;
    ret.cols = aa[0].length;

    foreach(i, row; aa)
    {
        enforce(row.length == ret.cols,
            "Cannot create a matrix from a jagged array of arrays.");

        foreach(j, elem; row) {
            ret[i, j] = elem;
        }
    }

    return ret;
}

unittest
{
    auto dense1 = matrix!real(3, 4);
    check (dense1.rows == 3  &&  dense1.cols == 4);
    check (isNaN(dense1[1,2]));

    auto dense2 = matrix(4, 3, 1.0);
    check (dense2.rows == 4  &&  dense2.cols == 3);
    check (dense2[1,2] == 1.0);

    auto upTri = matrix!(real, Storage.triangular)(3);
    check (upTri.rows == 3  &&  upTri.cols == 3);
    check (isNaN(upTri[0,2])  &&  upTri[2,0] == 0);

    auto loTri = matrix!(double, Storage.triangular, Triangle.lower)(3, 1.0);
    check (loTri.rows == 3  &&  loTri.cols == 3);
    check (loTri[0,2] == 0  &&  loTri[2,0] == 1.0);

}

/** Create a diagonal matrix from an array of elements.*/
MatrixView!(T, Storage.diagonal) diag(T)(T[] diagonalElements) pure nothrow
{
    return typeof(return)(
        diagonalElements,
        diagonalElements.length,
        diagonalElements.length
    );
}

/** A convenience function that creates a copy of the input matrix. */
MatrixView!(T, stor, tri) copy(T, Storage stor, Triangle tri)
    (const MatrixView!(T, stor, tri) m)
    pure
{
    MatrixView!(T, stor, tri) mcopy;
    mcopy.rows = m.rows;
    mcopy.cols = m.cols;
    mcopy.array = m.array.dup;
    return mcopy;
}

unittest
{
    auto a = matrix!double(2, 2);
    a[1,0] = 1.0;
    auto b = copy(a);
    b[1,0] = 2.0;
    check (b[1,0] == 2.0  &&  a[1,0] == 1.0);
}

/** This struct provides a matrix-like view of the contents of an
    array. In order to be compatible with LAPACK routines, it supports
    the following matrix representations (i.e. memory layouts).

    General_matrices:
    The elements of dense matrices are stored in column-major order.
    This means that if the wrapped array contains the elements
    ---
    a b c d e f g h i j k l
    ---
    then a 3x4 dense matrix view of this array looks like this:
    ---
    a d g j
    b e h k
    c f i l
    ---

    Triangular_matrices:
    Triangular matrices are required to be square. If the wrapped
    array contains the six elements
    ---
    a b c d e f
    ---
    then the resulting 3x3 upper and lower triangular matrix views
    will look like this, respectively:
    ---
    a b d         a 0 0
    0 c e   and   b d 0
    0 0 f         c e f
    ---

    Symmetric_matrices:
    Symmetric matrices are stored in the same way as triangular
    matrices. This means that for the array above, the corresponding
    symmetric matrix view will be
    ---
    a b d       a b c
    b c e   or  b d c
    d e f       c e f
    ---
    depending on whether the upper or lower triangle is stored.

    Hermitian_matrices:
    Hermitian matrices are not implemented yet.

    See_also:
    LAPACK User's Guide: Matrix storage schemes,
    $(LINK http://www.netlib.org/lapack/lug/node121.html)
*/
struct MatrixView (T, Storage stor = Storage.general,
    Triangle tri = Triangle.upper)
{
    enum Storage storage = stor;
    enum Triangle triangle = tri;

    /// false.
    enum isTransposed = false;

private:
    // Writing fully-qualified names in static ifs gets tiresome, so
    // we introduce a few flags.
    enum : bool
    {
        isGen   = (storage == Storage.general),
        isTri   = (storage == Storage.triangular),
        isUpTri = (isTri && triangle == Triangle.upper),
        isLoTri = (isTri && triangle == Triangle.lower),
        isSym   = (storage == Storage.symmetric),
        isUpSym = (isSym && triangle == Triangle.upper),
        isLoSym = (isSym && triangle == Triangle.lower),
        isDiag  = (storage == Storage.diagonal)
    }


    static if (!isGen)
    {
        static assert (isNumeric!T || isComplex!T,
            "MatrixView: Non-general matrices can only contain numeric values. "
           ~"Non-appropriate type given: "~T.stringof);
    }

    // The zero element in triangular and diagonal matrices.
    static if (isTri || isDiag)
    {
        T zero = Zero!T;
    }


public:
    /** The array that is wrapped by this MatrixView. */
    T[] array;


    /** The number of rows in the matrix. */
    size_t rows;

    /** The number of columns in the matrix. */
    size_t cols;

    /** The leading (row) dimension.  Included to support matrix slicing,
        currently just an alias to rows.
    */
    alias rows leading;



    /** Wrap a MatrixView with m rows around the given array.

        For general matrices, the number of columns in the matrix
        is set to a.length/m, whereas for triangular and symmetric
        matrices the number of columns is set equal to the number
        of rows.
    */
    this (T[] a, size_t m)  pure nothrow
    in
    {
        static if (isGen)  assert (a.length % m == 0);
    }
    body
    {
        static if (isGen)  this (a, m, a.length/m);
        else static if (isTri || isSym)  this(a, m, m);
    }



    /** Wrap a MatrixView with m rows and n columns around the given array.
        For a given set of a, m, and n, the following must be true for
        a general matrix:
        ---
        a.length >= m*n
        ---
        For triangular and symmetric matrices, there are two constraints:
        ---
        m == n
        a.length >= (n*n + n)/2
        ---
        These conditions are only checked in non-release builds.
    */
    this (T[] a, size_t m, size_t n) pure nothrow
    in
    {
        static if (isGen)
            assert (a.length >= m*n);
        else static if (isTri || isSym)
        {
            assert (m == n);
            assert (a.length >= (n*n + n)/2);
        }
    }
    body
    {
        array = a;
        rows = m;
        cols = n;
    }



    /** Return (a reference to) the element at row i, column j.

        Warning:
        For convenience, this method returns values by reference. This
        means that one can do stuff like this:
        ---
        m[1,2] += 3.14;
        ---
        Unfortunately, it also means that in a triangular or diagonal matrix one
        can change the zero element (which is common for all zero elements).
        ---
        assert ((m[1,0] == 0.0)  &&  m[2,0] == 0.0);
        m[1,0] += 3.14;
        assert (m[2,0] == 3.14);    // passes
        ---
        You are hereby warned.
    */
    ref T opIndex(size_t i, size_t j) pure nothrow
    in
    {
        assert (i < rows  &&  j < cols);
    }
    body
    {
        static if (isGen)
            return array.ptr[i + rows*j];
        else static if (isUpTri)
        {
            if (i <= j)  return array.ptr[i + (j*j+j)/2];
            else return zero;
        }
        else static if (isLoTri)
        {
            if (i >= j)  return array.ptr[i + ((rows+rows-j-1)*j)/2];
            else return zero;
        }
        else static if (isUpSym)
        {
            if (i <= j)  return array.ptr[i + (j*j+j)/2];
            else return array.ptr[j + (i*i+i)/2];
        }
        else static if (isLoSym)
        {
            if (i >= j)  return array.ptr[i + ((rows+rows-j-1)*j)/2];
            else return array.ptr[j + ((rows+rows-i-1)*i)/2];
        }
        else static if(isDiag)
        {
            if(i == j) return array[i];
            else return zero;
        }
        else static assert (false);
    }


    /** Assign a value to the element at row i, column j.

        Unlike opIndex(), this method checks that zero elements in
        a triangular or diagonal matrix aren't assigned to, but only in
        non-release builds.
    */
    T opIndexAssign(T value, size_t i, size_t j) nothrow
    in
    {
        assert (i < rows  &&  j < cols);
        static if (isUpTri)  assert (i <= j);
        static if (isLoTri)  assert (i >= j);
        static if (isDiag) assert(i == j);
    }
    body
    {
        static if (isGen)
            return array.ptr[i + rows*j] = value;
        else static if (isUpTri)
            return array.ptr[i + (j*j+j)/2] = value;
        else static if (isLoTri)
            return array.ptr[i + ((rows+rows-j-1)*j)/2] = value;
        else static if (isUpSym)
        {
            if (i <= j)  return array.ptr[i + (j*j+j)/2] = value;
            else  return array.ptr[j + (i*i+i)/2] = value;
        }
        else static if (isLoSym)
        {
            if (i >= j)  return array.ptr[i + ((rows+rows-j-1)*j)/2] = value;
            else return array.ptr[j + ((rows+rows-i-1)*i)/2] = value;
        }
        else static if(isDiag)
        {
            return array.ptr[i] = value;
        }
        else static assert (false);
    }

    ///
    Transposed!(typeof(this)) transpose() {
        return typeof(return)(this);
    }

    /**
    Prints a matrix in the default precision, which is 6 significant
    figures.
    */
    string toString()
    {
        return toString(6);
    }

    /**
    Prints a matrix using a user-specified number of significant figures.
    */
    string toString(uint precision)
    {
        auto app = appender!(string)();
        app.reserve(rows * cols * (precision + 6) + rows);
        auto formatStr = text('%', precision + 6, '.', precision, 'g');

        // The hard coded 6's are for the exponent, i.e. e-308, plus a space.
        foreach(i; 0..rows)
        {
            foreach(j; 0..cols)
            {
                app.put(
                    format(formatStr, this[i, j])
                );
            }

            app.put('\n');
        }

        return app.data;
    }

    /**
    Generates an expression template for adding or subtracting rhs to/from this.
    rhs must have the same dimensions, but need not have the same storage.

    For details on expression templates, see scid.matrixops.AddSubExpr.
    */
    auto opBinary(string op, M)(M rhs)
    if((isMatrixView!M || isTransposed!M) && (op == "+" || op == "-"))
    {
        return binaryImpl!op(this, rhs);
    }

    /**
    Test two matrices for equality.  They need not have the exact same type
    or same storage.
    */
    bool opEquals(M)(const M rhs) const
    if(isMatrixView!M)
    {
        static if(rhs.storage == this.storage)
        {
            return rows == rhs.rows && cols == rhs.cols && array == rhs.array;
        }
        else
        {
            static assert(0,
                "Equality testing with different storage not implemented yet.");
        }
    }

    /// ditto
    bool opEquals(Templ)(Templ exprTempl) const
    if(is(typeof(exprTempl.evaluate())))
    {
        return opEquals(exprTempl.evaluate());
    }
}

unittest
{
    alias MatrixView!real GeneralMatrix;
    real[] g = [1.0L, 2, 3, 4, 5, 6];

    auto gm1 = GeneralMatrix(g, 2);
    auto gm2 = GeneralMatrix(g, 3);

    check (gm1.cols == 3);
    check (gm2.cols == 2);

    check (gm1[1,0] == 2);
    check (gm2[1,0] == 2);

    check (gm1[1,1] == 4);
    check (gm2[1,1] == 5);

    gm2[1,1] += 1; check (gm2[1,1] == 6);
    gm2[1,1] = 10; check (gm2[1,1] == 10);


    alias MatrixView!(real, Storage.triangular) UTMatrix;
    real[] u = [1.0, 2, 3, 4, 5, 6];

    auto um1 = UTMatrix(u, 3);
    check (um1.cols == 3);
    check (um1[1,0] == 0.0);
    check (um1[1,1] == 3.0);
    um1[0,2] += 3; check (u[3] == 7);
    um1[2,2] = 10; check (u[5] == 10);


    alias MatrixView!(real, Storage.triangular, Triangle.lower) LTMatrix;
    real[] l = [1.0, 2, 3, 4, 5, 6];

    auto lm1 = LTMatrix(l, 3);
    check (lm1.cols == 3);
    check (lm1[0,1] == 0.0);
    check (lm1[1,1] == 4.0);
    lm1[2,0] += 4; check (l[2] == 7);
    lm1[2,2] = 10; check (l[5] == 10);


    alias MatrixView!(real, Storage.symmetric) USMatrix;
    real[] us = [1.0, 2, 3, 4, 5, 6];

    auto usm1 = USMatrix(us, 3);
    check (usm1.cols == 3);
    check (usm1[1,2] == 5.0);
    foreach (i; 0 .. usm1.rows)
        foreach (j; 0 .. i)
            check (usm1[i,j] == usm1[j,i]);
    usm1[0,2] += 3; check (usm1[2,0] == 7);
    usm1[1,2] = 10; check (usm1[2,1] == 10);


    alias MatrixView!(real, Storage.symmetric, Triangle.lower) LSMatrix;
    real[] ls = [1.0, 2, 3, 4, 5, 6];

    auto lsm1 = LSMatrix(ls, 3);
    check (lsm1.cols == 3);
    check (lsm1[1,2] == 5.0);
    foreach (i; 0 .. lsm1.rows)
        foreach (j; 0 .. i)
            check (lsm1[i,j] == lsm1[j,i]);
    lsm1[0,2] += 3; check (lsm1[2,0] == 6);
    lsm1[1,2] = 10; check (lsm1[2,1] == 10);

    real[] diags = [1.0, 2, 3, 4, 5, 6];
    auto dmat = diag(diags);
    check(dmat.rows == 6);
    check(dmat.cols == 6);
    check(dmat[0, 2] == 0);
    check(dmat[1, 1] == 2);
}

/** Evaluates to true if the given type is an instantiation of
    MatrixView. Optionally test the element type and/or storage
    scheme.
*/
template isMatrixView(MatrixT)
{
    static if (is(MatrixT E == MatrixView!(E, S, T), int S, char T))
    {
        enum bool isMatrixView = true;
    }
    else enum bool isMatrixView = false;
}

version(unittest)
{
    static assert (isMatrixView!(MatrixView!int));
    static assert (!isMatrixView!int);
}


/// ditto
template isMatrixView(MatrixT, ElemT)
{
    static if (is(MatrixT E == MatrixView!(E, S, T), int S, char T))
    {
        enum bool isMatrixView = is(E == ElemT);
    }
    else enum bool isMatrixView = false;
}

version(unittest)
{
    static assert (isMatrixView!(MatrixView!int, int));
    static assert (!isMatrixView!(MatrixView!int, float));
}


/// ditto
template isMatrixView(MatrixT, Storage stor)
{
    static if (is(MatrixT E == MatrixView!(E, S, T), int S, char T))
    {
        enum bool isMatrixView = S == stor;
    }
    else enum bool isMatrixView = false;
}

version(unittest)
{
    static assert (isMatrixView!(
        MatrixView!(int, Storage.triangular),
        Storage.triangular));
    static assert (!isMatrixView!(
        MatrixView!(int, Storage.triangular),
        Storage.symmetric));
}


/// ditto
template isMatrixView(MatrixT, ElemT, Storage stor)
{
    static if (is(MatrixT E == MatrixView!(E, S, T), int S, char T))
    {
        enum bool isMatrixView = is(E == ElemT)  &&  S == stor;
    }
    else enum bool isMatrixView = false;
}

version(unittest)
{
    static assert (isMatrixView!(
        MatrixView!(int, Storage.triangular),
        int, Storage.triangular));
    static assert (!isMatrixView!(
        MatrixView!(int, Storage.triangular),
        float, Storage.triangular));
    static assert (!isMatrixView!(
        MatrixView!(int, Storage.triangular),
        int, Storage.symmetric));
}

template MatrixType(M)
if(isMatrixView!M)
{
    alias typeof(M.init.array[0]) MatrixType;
}

unittest
{
    static assert(is(MatrixType!(MatrixView!double) == double));
}

/**
Holds a transposed matrix, i.e. the columns become the rows and the rows
become the columns.
*/
struct Transposed(Matrix)
if(isMatrixView!Matrix)
{
    // For convenience/typing saving
    enum Storage storage = matrix.storage;
    enum Triangle triangle = matrix.triangle;

    /// The type of the contents of Matrix.
    alias typeof(matrix[0, 0]) E;

    /// true
    enum bool isTransposed = true;

    /// The underlying matrix.
    Matrix matrix;

    /// Returns matrix.cols
    size_t rows() @property
    {
        return matrix.cols;
    }

    /// Returns matrix.rows
    size_t cols() @property
    {
        return matrix.rows;
    }

    ///
    ref E opIndex(size_t i, size_t j)
    {
        return matrix[j, i];
    }

    ///
    E opIndexAssign(E val, size_t i, size_t j)
    {
        return matrix[j, i] = E;
    }

    ///
    auto opBinary(string op, M)(M rhs)
    if((isMatrixView!M || isTransposed!M) && (op == "+" || op == "-"))
    {
        return binaryImpl!op(this, rhs);
    }

    ///
    Matrix transpose() {
        return matrix;
    }
}

private auto binaryImpl(string op, M1, M2)(M1 lhs, M2 rhs)
if(op == "+" || op == "-")
{
    auto lhsWithSign = plusMinusMatrix!('+')(lhs);
    auto rhsWithSign = plusMinusMatrix!(op[0])(rhs);

    AddSubExpr!(typeof(lhsWithSign), typeof(rhsWithSign)) ret;
    ret.matrices[0] = lhsWithSign;
    ret.matrices[1] = rhsWithSign;
    return ret;
}

template Transposed(T)
if(isTransposed!T)
{
    alias typeof(typeof(T.init.matrix)) Transposed;
}

template isTransposed(M)
{
    // Just look at the duck interface.
    enum bool isTransposed = is(typeof(M.isTransposed)) && M.isTransposed;
}

template CommonMatrix(M...)
{
    alias CommonMatrixImpl!(M).ret CommonMatrix;
}

private template CommonMatrixImpl(M...)
{
    alias CommonType!(
        typeof(M[0].init.array[0]),
        typeof(M[1].init.array[0])
    ) E;

    static if(M.length == 1)
    {
        alias M ret;
    }
    else static if(M.length == 2)
    {

        static if(M[0].storage == Storage.diagonal && M[1].storage ==
        Storage.diagonal)
        {
            alias MatrixView!(E, Storage.diagonal) ret;
        }
        else static if(M[0].storage == Storage.triangular
        && M[1].storage == Storage.diagonal)
        {
            alias MatrixView!(E, Storage.triangular, M[0].triangle) ret;
        }
        else static if(M[1].storage == Storage.triangular
        && M[0].storage == Storage.diagonal)
        {
            alias MatrixView!(E, Storage.triangular, M[1].triangle) ret;
        }
        else static if(M[0].storage == Storage.triangular
        && M[1].storage == Storage.triangular && M[0].triangle == M[1].triangle)
        {
            alias MatrixView!(E, Storage.triangular, M[1].triangle) ret;
        }
        else
        {
            alias MatrixView!(E, Storage.general) ret;
        }
    }
    else
    {
        alias CommonMatrix!(
            CommonMatrix!(M[0], M[1]), M[2..$]
        ) ret;
    }
}

unittest
{
    static assert(is(CommonMatrix!(
        MatrixView!(double, Storage.diagonal),
        MatrixView!(float, Storage.diagonal)) ==
        MatrixView!(double, Storage.diagonal)
    ));

    static assert(is(CommonMatrix!(
        MatrixView!(ubyte, Storage.general),
        MatrixView!(float, Storage.diagonal)) ==
        MatrixView!(float, Storage.general)
    ));

    static assert(is(CommonMatrix!(
        MatrixView!(double, Storage.triangular, Triangle.upper),
        MatrixView!(float, Storage.triangular, Triangle.lower)) ==
        MatrixView!(double, Storage.general)
    ));

    static assert(is(CommonMatrix!(
        MatrixView!(double, Storage.triangular, Triangle.upper),
        MatrixView!(float, Storage.triangular, Triangle.upper)) ==
        MatrixView!(double, Storage.triangular, Triangle.upper)
    ));
}

// CTFE function.
bool allSameStorage(M...)()
{
    if(M.length < 2) return true;

    foreach(i, m; M[0..$ - 1])
    {
        static if(M[0].storage != M[1].storage)
        {
            return false;
        }
        else static if(M[0].storage != Storage.triangular)
        {
            continue;
        }
        else static if(M[0].triangle == M[1].triangle)
        {
            continue;
        }
        else
        {
            return false;
        }
    }

    return true;
}

// CTFE function
bool noneTransposed(M...)()
{
    foreach(m; M)
    {
        if(m.isTransposed) return false;
    }

    return true;
}
