module scid.matrixops;

import std.typetuple;
import scid.matrix;

struct PlusMinusMatrix(T, Storage stor, Triangle tri, char sign_)
{
    MatrixView!(T, stor, tri) matrix;
    alias matrix this;
    enum char sign = sign_;

    static if(sign == '-')
    ref T opIndex(size_t i, size_t j)
    {
        return -matrix[i, j];
    }
}

template isPlusMinusMatrix(T)
{
    // Look at duck interface.
    enum bool isPlusMinusMatrix =
        is(typeof(T.sign == char)) &&
        is(typeof(T.matrix.storage == Storage)) &&
        is(typeof(T.matrix.triangle == Triangle)) &&
        isArray!(typeof(T.matrix.array));
}

template BaseMatrix(T)
{
    alias typeof(T.init.matrix) BaseMatrix;
}

struct AddSubTempl(M...)
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
    {
        static if(allSameStorage!(M))
        {
            result.rows = matrices[0].rows;
            result.cols = matrices[0].cols;
            result.array = newVoid!(typeof(result.array[0]))
                (matrices[0].array.length);

            static if(M.length == 2 &&
            is(typeof(M[0].init.array) == typeof(M[1].init.array))
            {
                // Use array ops.
                mixin("result.array[] = " ~ matrices[0].sign ~
                    "matrices[0].array " ~ matrices[1].sign ~
                    "matrices[1].array;");
            }
            else
            {
                // Use fused manual loop.
                foreach(i; 0..result.array.length)
                {
                    mixin("results.array[i] = " ~
                        makeAddSubExpr(matrices) ~ ';');
                }
            }
        }
        else
        {
            static assert(0, "Matrix addition/subtraction with heterogeneous "
                ~ " storage not implemented yet.");
        }

        return result;
    }

    alias evaluate this;
}

import std.conv;

// CTFE function to generate the code for an iteration of the addition/
// subtraction loop.
private string makeAddSubExpr(M...)(M matrices)
{
    string ret;
    foreach(i, m; matrices)
    {
        ret ~= ' ' ~ sign(m);
        ret ~= ' ' ~ "matrices[" ~ to!string(i) ~ "][i]";
    }

    return ret;
}
