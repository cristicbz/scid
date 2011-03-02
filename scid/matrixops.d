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
    {
        static if(allSameStorage!(M)())
        {
            result.rows = matrices[0].rows;
            result.cols = matrices[0].cols;
            result.array = newVoid!(typeof(result.array[0]))
                (matrices[0].array.length);

            static if(sameArrayTypes!M())
            {
                // Use array ops.
                enum toMixIn = makeArrayAddSubExpr!M();
                result.array[] = mixin(toMixIn);
            }
            else
            {
                // Use fused manual loop.
                foreach(i; 0..result.array.length)
                {
                    mixin("result.array[i] = " ~
                        makeAddSubExpr!M() ~ ';');
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
    if(isMatrixView!M2 && (op == "+" || op == "-"))
    {
        auto rhsWithSign = plusMinusMatrix!(op[0])(rhs);
        return AddSubExpr!(M, typeof(rhsWithSign))(matrices, rhsWithSign);
    }

    auto opBinaryRight(string op, M2)(M2 lhs)
    if(isMatrixView!M2 && (op == "+" || op == "-"))
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

// CTFE function to generate the code for an iteration of the addition/
// subtraction loop.
private string makeAddSubExpr(M...)()
{
    string ret;
    foreach(i, m; M)
    {
        if(i > 0 || m.sign == '-')
        {
            ret ~= " " ~ m.sign;
        }

        ret ~= ' ' ~ "matrices[" ~ to!string(i) ~ "].array[i]";
    }

    return ret;
}

// CTFE function to generate the code for an iteration of the addition/
// subtraction array ops.
private string makeArrayAddSubExpr(M...)()
{
    string ret;
    foreach(i, m; M)
    {
        if(i > 0 || m.sign == '-')
        {
            ret ~= " " ~ m.sign;
        }

        ret ~= ' ' ~ "matrices[" ~ to!string(i) ~ "].array[]";
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
