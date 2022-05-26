#include "boost.h"

Matrix transpose(Matrix& A)
{
    Matrix B(A.size2(), A.size1());
    for (size_t i = 0; i < A.size1(); i++)
        for (size_t j = 0; j < A.size2(); j++)
            B(j, i) = A(i, j);
    return B;
}

triLowerMatrix choletsky(symLowerMatrix& A)
{
    int n = A.size1();
    double S1, S2;
    triLowerMatrix L(n, n);

    L(0, 0) = std::sqrt(A(0,0));

    for (int i = 1; i < n; i++)
        L(i,0) = A(i,0) / L(0, 0);

    for (int k = 1; k < n; k++) {
        S1 = 0;

        for (int i = 0; i < k; i++)
            S1 += L(k,i) * L(k, i);

        L(k, k) = std::sqrt(A(k,k) - S1);

        for (int j = k + 1; j < n; j++) {
            S2 = 0;

            for (int i = 0; i < k; i++)
                S2 += L(k, i) * L(j, i);

            L(j,k) = (A(j,k) - S2) / L(k,k);
        }
    }
    return L;
}

triUpperMatrix transpose(triLowerMatrix& A)
{
    triUpperMatrix B(A.size1(), A.size1());
    for (size_t i = 0; i < A.size1(); i++)
        for (size_t j = i; j < A.size1(); j++)
            B(i, j) = A(j, i);
    return B;
}

symLowerMatrix triProd(triUpperMatrix& A, triLowerMatrix& B)
{
    size_t n = A.size1();
    symLowerMatrix C(n, n);
    for (size_t j = 0; j < n; j++)
        for (size_t i = j; i < n; i++)
        {
            double val = 0;
            for (size_t k = i; k < n; k++)
                val += A(i, k) * B(k, j);
            C(i, j) = val;
        }
    return C;
}
