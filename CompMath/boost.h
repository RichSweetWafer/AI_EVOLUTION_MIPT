#ifndef BOOST_H
#define BOOST_H
// Random number generator
#include <boost/random.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/constants/constants.hpp>

typedef boost::mt19937 RNGType;
#include <QDebug>
// Matrices
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

typedef boost::numeric::ublas::lower lower;
typedef boost::numeric::ublas::upper upper;
typedef boost::numeric::ublas::matrix<double> Matrix;
typedef boost::numeric::ublas::symmetric_matrix<double, lower> symLowerMatrix;
typedef boost::numeric::ublas::symmetric_matrix<double, upper> symUpperMatrix;
typedef boost::numeric::ublas::triangular_matrix<double, lower> triLowerMatrix;
typedef boost::numeric::ublas::triangular_matrix<double, upper> triUpperMatrix;


// Matrix operations
using boost::numeric::ublas::trans;
using boost::numeric::ublas::prod;

Matrix transpose(Matrix& A);

triLowerMatrix choletsky(symLowerMatrix& A);
triUpperMatrix transpose(triLowerMatrix& A);
symLowerMatrix triProd(triUpperMatrix& A, triLowerMatrix& B);
triLowerMatrix invTriLower(triLowerMatrix& A);
symLowerMatrix inverse(symLowerMatrix& A);
//symLowerMatrix inverseL(triLowerMatrix& L);
#endif // BOOST_H
