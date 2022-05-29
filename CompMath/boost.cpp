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
    int n = A.size1();
    symLowerMatrix C(n, n);
    for (int j = 0; j < n; j++)
        for (int i = j; i < n; i++)
        {
            double val = 0;
            for (int k = i; k < n; k++)
                val += A(i, k) * B(k, j);
            C(i, j) = val;
        }
    return C;
}

triLowerMatrix invTriLower(triLowerMatrix& A)
{
    int size = A.size1();
    triLowerMatrix B(size, size);
    for (int i = 0; i < size; i++)
        B(i, i) = 1 / A(i, i);
    for (int i = 1; i < size; i++)
        for (int j = 0; j < i; j++)
        {
            B(i, j) = - double(1 / A(i, i));
            double sum = 0;
            for (int k = 0; k < i; k++)
            {
                if (k < j || k > i)
                    continue;
                sum += B(k, j)*A(i, k);
            }
            B(i, j) *= sum;
        }
    return B;
}

symLowerMatrix inverse(symLowerMatrix& A)
{
    triLowerMatrix tri = choletsky(A);
    triLowerMatrix invtri = invTriLower(tri);
    triUpperMatrix invtriTr = transpose(invtri);
    return triProd(invtriTr, invtri);
}

//symLowerMatrix inverseL(triLowerMatrix& L)
//{
//    int size = L.size1();
//    triLowerMatrix B(size, size);
//    for (int i = 0; i < size; i++)
//        B(i, i) = 1 / L(i, i);
//    for (int i = 1; i < size; i++)
//        for (int j = 0; j < i; j++)
//        {
//            B(i, j) = - double(1 / L(i, i));
//            double sum = 0;
//            for (int k = 0; k < i; k++)
//            {
//                if (k < j || k > i)
//                    continue;
//                sum += B(k, j)*L(i, k);
//            }
//            B(i, j) *= sum;
//        }
//    triUpperMatrix C(B.size1(), B.size1());
//    for (size_t i = 0; i < B.size1(); i++)
//        for (size_t j = i; j < B.size1(); j++)
//            C(i, j) = B(j, i);
//    return triProd(C, B);
//    // for invL
//    triUpperMatrix B(L.size1(), L.size1());
//    for (size_t i = 0; i < L.size1(); i++)
//        for (size_t j = i; j < L.size1(); j++)
//            B(i, j) = L(j, i);

//    int n = B.size1();
//    symLowerMatrix C(n, n);
//    for (int j = 0; j < n; j++)
//        for (int i = j; i < n; i++)
//        {
//            double val = 0;
//            for (int k = i; k < n; k++)
//                val += B(i, k) * L(k, j);
//            C(i, j) = val;
//        }
//    return C;
//}
