#ifndef CPT_LINALG_HPP
#define CPT_LINALG_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>

#include "matrix.hpp"

namespace cpt
{

    extern Matrix<double,2> transpose(const Matrix<double,2>&);

    extern void solve_Gauss_Jordan(Matrix<double,2>& A, Matrix<double,2>& B);

    extern void solve_LU_decompose(Matrix<double,2>& A, Matrix<double,2>& B);

    extern Matrix<double,1> solve_eigen_symmetric(Matrix<double,2>& A);

    extern Matrix<double,1> solve_eigen_generalized(Matrix<double,2>& A,
                                                    Matrix<double,2>& S);

    template<class T> Matrix<T,2> operator+(
        const Matrix<T,2>& m1, const Matrix<T,2>&m2)
    {
        int rows = m1.dim1();
        int cols = m1.dim2();
        if (rows != m2.dim1() || cols != m2.dim2())
            error("linalg: matrix dimension mismatch");
        Matrix<T,2> m1_plus_m2(rows, cols);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                m1_plus_m2[i][j] = m1[i][j] + m2[i][j];
        return m1_plus_m2;
    }

    template<class T> Matrix<T,2> operator-(
        const Matrix<T,2>& m1, const Matrix<T,2>&m2)
    {
        int rows = m1.dim1();
        int cols = m1.dim2();
        if (rows != m2.dim1() || cols != m2.dim2())
            error("linalg: matrix dimension mismatch");
        Matrix<T,2> m1_minus_m2(rows, cols);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                m1_minus_m2[i][j] = m1[i][j] - m2[i][j];
        return m1_minus_m2;
    }

    template<class T> Matrix<T,2> operator*(
        const Matrix<T,2>& m1, const Matrix<T,2>&m2)
    {
        int n = m1.dim2();
        if (n != m2.dim1())
            error("linalg: matrix dimension mismatch");
        int rows = m1.dim1();
        int cols = m2.dim2();
        Matrix<T,2> product(rows, cols);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                for (int k = 0; k < n; k++)
                    product[i][j] += m1[i][k] * m2[k][j];
        return product;
    }

    template<class T> Matrix<T,1> operator+(
        const Matrix<T,1>& v1, const Matrix<T,1>& v2)
    {
        if (v1.dim1() != v2.dim1())
            error("linalg: matrix dimension mismatch");
        Matrix<T,1> v1_plus_v2(v1.dim1());
        for (int i = 0; i < v1.dim1(); i++)
            v1_plus_v2[i] = v1[i] + v2[i];
        return v1_plus_v2;
    }

    template<class T> Matrix<T,1> operator-(
        const Matrix<T,1>& v1, const Matrix<T,1>& v2)
    {
        if (v1.dim1() != v2.dim1())
            error("linalg: matrix dimension mismatch");
        Matrix<T,1> v1_minus_v2(v1.dim1());
        for (int i = 0; i < v1.dim1(); i++)
            v1_minus_v2[i] = v1[i] - v2[i];
        return v1_minus_v2;
    }

    template<class T> T dot(const Matrix<T,1>& v1, const Matrix<T,1>& v2)
    {
        T dot = 0;
        for (int i = 0; i < v1.size(); i++)
            dot += v1[i] * v2[i];
        return dot;
    }

    // FFT Routines

    extern void fft(Matrix<std::complex<double>,1>&,
                    const bool inverse=false);
    inline void fft_inv(Matrix<std::complex<double>,1>& v) { fft(v, true); }

    extern void fft(Matrix<std::complex<double>,2>&,
                    const bool inverse=false);
    inline void fft_inv(Matrix<std::complex<double>,2>& m) { fft(m, true); }

    // Least Squares and Chi-Square Fits

    extern void least_squares_fit(      // makes a linear least-squares fit
        const Matrix<double,1>& x,      // vector of x values - input
        const Matrix<double,1>& y,      // vector of y values - input
        double& a,                      // fitted intercept - output
        double& b,                      // fitted slope - output
        double& sigma_a,                // estimated error in intercept - output
        double& sigma_b,                // estimated error in slope - output
        double& sigma);                 // estimated error bar in y

    extern void chi_square_fit(         // makes a linear chi-square fit
        const Matrix<double,1>& x,      // vector of x values - input
        const Matrix<double,1>& y,      // vector of y values - input
        const Matrix<double,1>& err,    // vector of y error values - input
        double& a,                      // fitted intercept - output
        double& b,                      // fitted slope - output
        double& sigma_a,                // estimated error in intercept - output
        double& sigma_b,                // estimated error in slope - output
        double& chi_square);            // minimized chi-square sum - output

} /* namespace cpt */

#endif /* CPT_LINALG_HPP  */
