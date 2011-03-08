module scid.matrixops;

import std.conv;
import std.traits;
import std.typetuple;
import scid.matrix;
import scid.core.memory;

version(unittest)
{
    import std.stdio;
}

package struct PlusMinusMatrix(T, Storage stor, Triangle tri, char sign_)
{
    MatrixView!(T, stor, tri) matrix;
    alias matrix this;
    enum char sign = sign_;

    static if(sign == '-')
    T opIndex(size_t i, size_t j)
    {
        return -matrix[i, j];
    }
}

package auto plusMinusMatrix(char sign, T, Storage stor, Triangle tri)
(MatrixView!(T, stor, tri) mat)
{
    return PlusMinusMatrix!(T, stor, tri, sign)(mat);
}

template isPlusMinusMatrix(T)
{
    // Look at duck interface.
    enum bool isPlusMinusMatrix =
        is(typeof(T.sign) == char) &&
        is(typeof(T.matrix.storage == Storage)) &&
        is(typeof(T.matrix.triangle == Triangle)) &&
        isArray!(typeof(T.matrix.array));
}

template BaseMatrix(T)
{
    alias typeof(T.init.matrix) BaseMatrix;
}

// CTFE function.
bool sameArrayTypes(M...)()
{
    foreach(i, m; M[0..$ - 1])
    {
        if(!is(typeof(m.init.array) == typeof(M[i + 1].init.array)))
        {
            return false;
        }
    }

    return true;
}

private template Negate(M)
{
    static if(M.sign == '-')
    {
        alias PlusMinusMatrix!(typeof(M.init.array[0]),
            M.storage, M.triangle, '+') Negate;
    }
    else
    {
        static assert(M.sign == '+');
        alias PlusMinusMatrix!(typeof(M.init.array[0]),
            M.storage, M.triangle, '-') Negate;
    }
}

struct AddSubExpr(M...)
if(allSatisfy!(isPlusMinusMatrix, M))
{
    M matrices;

    alias CommonMatrix!(staticMap!(BaseMatrix, M)) ResultType;
    ResultType result;

    ResultType evaluate()
    {
        // The point of this coding style is to get the common case to
        // be inlined.
        return (result.rows > 0) ? result : evaluateImpl();
    }

    ResultType evaluateImpl()
    in
    {
        foreach(m; matrices[1..$])
        {
            assert(m.rows == matrices[0].rows,
                "Matrices to be added/subtracted must have same number of rows.");
            assert(m.cols == matrices[0].cols,
                "Matrices to be added/subtracted must have same number of cols.");
        }
    }
    body
    {
        result.rows = matrices[0].rows;
        result.cols = matrices[0].cols;
        alias typeof(result.array[0]) E;

        static if(result.storage == Storage.general)
        {
            result.array = newVoid!(E)(matrices[0].array.length);
        }
        else static if(result.storage == Storage.symmetric ||
        result.storage == Storage.triangular)
        {
            // Symmetric and triangular matrices have to be square.
            size_t n = result.rows;
            assert(result.cols == n);
            result.array = newVoid!(E)((n + n * n) / 2);
        }
        else
        {
            // Diagonal matrices have to be square, too.
            static assert(result.storage == Storage.diagonal);
            size_t n = result.rows;
            assert(result.cols == n);
            result.array = newVoid!(E)(n);
        }

        static if(allSameStorage!(M)() && noneTransposed!(M)())
        {
            static if(sameArrayTypes!M())
            {
                // Use array ops.
                enum toMixIn = makeAddSubExpr!(everything, M)(".array[]");
                result.array[] = mixin(toMixIn);
            }
            else
            {
                // Use fused manual loop.
                foreach(i; 0..result.array.length)
                {
                    result.array[i] = mixin(
                        makeAddSubExpr!(everything, M)(".array[i]")
                    );
                }
            }
        }
        else
        {
            if(isSquare(result))
            {
                // Then we can have mixed storage types and things get delicate.
                // Add the diagonal elements first.
                size_t n = result.rows;
                foreach(i; 0..n)
                {
                    result[i, i] = mixin(
                        makeAddSubExpr!(everything, M)("[i, i]")
                    );
                }

                // Now do the lower triangle if necessary.  Handle the
                // symmetric case here arbitrarily instead of in the upper
                // triangle.
                static if(result.storage == Storage.general ||
                (result.storage == Storage.triangular &&
                 result.triangle == Triangle.lower) ||
                 result.storage == Storage.symmetric)
                 {
                     foreach(i; 0..n) foreach(j; 0..i)
                     {
                         result[i, j] = mixin(
                            makeAddSubExpr!(hasLowerTriangle, M)("[i, j]")
                        );
                     }
                 }

                 // Now do the upper triangle if necessary.  In the symmetric
                 // storage case we're done already.
                static if(result.storage == Storage.general ||
                (result.storage == Storage.triangular &&
                 result.triangle == Triangle.upper))
                 {
                     foreach(i; 0..n) foreach(j; i)
                     {
                         result[i, j] = mixin(
                            makeAddSubExpr!(hasUpperTriangle, M)("[i, j]")
                        );
                     }
                 }
            }
            else
            {
                // Then all the storage types should be general.  Double check.
                static assert(result.storage == Storage.general);
                foreach(Mat; M) static assert(Mat.storage == Storage.general);

                // We want to iterate with the grain for as many matrices
                // as possible.  Count how many matrices are transposed vs.
                // not.
                enum nTransposed = .nTransposed!(M, ResultType)();
                static if(nTransposed > (M.length + 1) / 2)
                {
                    // Then go rows first.
                    foreach(i; 0..result.rows) foreach(j; 0..result.cols)
                    {
                        result[i, j] = mixin(
                            makeAddSubExpr!(dummyFilter, M)("[i, j]")
                        );
                    }
                }
                else
                {
                    // Then go columns first.
                    foreach(j; 0..result.cols) foreach(i; 0..result.rows)
                    {
                        result[i, j] = mixin(
                            makeAddSubExpr!(dummyFilter, M)("[i, j]")
                        );
                    }
                }
            }
        }

        return result;
    }

    // Returns the negative of this expression template.
    auto negate()
    {
        AddSubExpr!(staticMap!(Negate, M)) ret;
        foreach(i, m; M)
        {
            ret.matrices[i].matrix = this.matrices[i].matrix;
        }

        return ret;
    }

    alias evaluate this;

    auto opBinary(string op, M2)(M2 rhs)
    if((isMatrixView!M2 || isTransposed!M2)  && (op == "+" || op == "-"))
    {
        auto rhsWithSign = plusMinusMatrix!(op[0])(rhs);
        return AddSubExpr!(M, typeof(rhsWithSign))(matrices, rhsWithSign);
    }

    auto opBinaryRight(string op, M2)(M2 lhs)
    if((isMatrixView!M2 || isTransposed!M2) && (op == "+" || op == "-"))
    {
        auto lhsWithSign = plusMinusMatrix!('+')(lhs);
        alias typeof(negate().matrices) NegateType;

        static if(op == "+")
        {
            return AddSubExpr!(typeof(lhsWithSign), M)
                (lhsWithSign, matrices);
        }
        else
        {
            return AddSubExpr!(typeof(lhsWithSign), NegateType)
                (lhsWithSign, negate().matrices);
        }
    }
}

// These templates are filters for various matrix storage types.
private template everything(T)
{
    enum bool everything = true;
}

private template hasLowerTriangle(T)
{
    enum bool hasLowerTriangle = T.storage == Storage.symmetric ||
        (T.storage == storage.triangular && T.triangle == Triangle.lower) ||
        T.storage == Storage.general;
}

private template hasUpperTriangle(T)
{
    enum bool hasUpperTriangle = T.storage == Storage.symmetric ||
        (T.storage == storage.triangular && T.triangle == Triangle.upper) ||
        T.storage == Storage.general;
}

// Same, but only does matrices that pass the filter.
private string makeAddSubExpr(alias filterFun, M...)(string op)
{
    string ret;
    foreach(i, m; M)
    {
        if(!filterFun!m) continue;

        if(i > 0 || m.sign == '-')
        {
            ret ~= " " ~ m.sign;
        }

        ret ~= ' ' ~ "matrices[" ~ to!string(i) ~ "]" ~ op;
    }

    return ret;
}

unittest {
    auto mat1 = matrix([[2.0, 1.0], [3.0, 4]]);
    auto mat2 = matrix([[8., 6], [7.0, 5]]);
    auto res1 = mat1 + mat2;

    assert(res1 == matrix([[10, 7], [10, 9]]));
    auto mat3 = matrix([[1, 2], [3, 4]]);
    auto res2 = res1 + mat3;
    writeln(res2);
    writeln(mat1 + mat2 + mat3);

    writeln(mat3 - res1);
}

void main() {}
