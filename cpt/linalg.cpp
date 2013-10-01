#include "linalg.hpp"
#include <cassert>

using namespace std;

namespace cpt {

Matrix<double,2> transpose(
    const Matrix<double,2>& A)
{
    Matrix<double,2> A_transpose(A.dim2(), A.dim1());
    for (int i = 0; i < A.dim1(); i++)
        for (int j = 0; j < A.dim2(); j++)
            A_transpose[j][i] = A[i][j];
    return A_transpose;
}

double inverse(Matrix<double,2>  A, Matrix<double,2>& Ainv) 
// Input
//    A    -    Matrix A (N by N)
// Outputs
//   Ainv  -    Inverse of matrix A (N by N)
//  determ -    Determinant of matrix A (return value)
{

  int N = A.dim1();
  assert( N == A.dim2() );
  
  Ainv = A;  // Copy matrix to ensure Ainv is same size
    
  int i, j, k;
  Matrix<double,1> scale(N);
  Matrix<double,2> b(N,N);       // Scale factor and work array
  int *index;  index = new int [N];

  //* Matrix b is initialized to the identity matrix
  b = 0.0;
  for( i=0; i<N; i++ )
    b(i,i) = 1.0;

  //* Set scale factor, scale(i) = max( |a(i,j)| ), for each row
  for( i=0; i<N; i++ ) {
    index[i] = i;                         // Initialize row index list
    double scalemax = 0.;
    for( j=0; j<N; j++ ) 
      scalemax = (scalemax > fabs(A(i,j))) ? scalemax : fabs(A(i,j));
    scale(i) = scalemax;
  }  

  //* Loop over rows k = 0, ..., (N-2)
  int signDet = 1;
  for( k=0; k<N-1; k++ ) {
        //* Select pivot row from max( |a(j,k)/s(j)| )
    double ratiomax = 0.0;
        int jPivot = k;
    for( i=k; i<N; i++ ) {
      double ratio = fabs(A(index[i],k))/scale(index[i]);
      if( ratio > ratiomax ) {
        jPivot=i;
        ratiomax = ratio;
      }
    }
        //* Perform pivoting using row index list
        int indexJ = index[k];
        if( jPivot != k ) {               // Pivot
      indexJ = index[jPivot];
      index[jPivot] = index[k];   // Swap index jPivot and k
      index[k] = indexJ;
          signDet *= -1;                          // Flip sign of determinant
        }
        //* Perform forward elimination
    for( i=k+1; i<N; i++ ) {
      double coeff = A(index[i],k)/A(indexJ,k);
      for( j=k+1; j<N; j++ )
        A(index[i],j) -= coeff*A(indexJ,j);
      A(index[i],k) = coeff;
      for( j=0; j<N; j++ ) 
        b(index[i],j) -= A(index[i],k)*b(indexJ,j);
    }
  }
  //* Compute determinant as product of diagonal elements
  double determ = signDet;         // Sign of determinant
  for( i=0; i<N; i++ )
        determ *= A(index[i],i);
  //* Perform backsubstitution
  for( k=0; k<N; k++ ) {
    Ainv(N-1,k) = b(index[N-1],k)/A(index[N-1],N-1);
    for( i=N-1; i>=0; i--) {
      double sum = b(index[i],k);
      for( j=i+1; j<N; j++ )
        sum -= A(index[i],j)*Ainv(j,k);
      Ainv(i,k) = sum/A(index[i],i);
    }
  }
  delete [] index;      // Release allocated memory
  return( determ );        
}

void solve_Gauss_Jordan(Matrix<double,2>& A, Matrix<double,2>& B) {

    int n = A.dim1();
    if (A.dim2() != n)
        error("linalg: matrix dimension mismatch");
    if (B.dim1() != n)
        error("linalg: matrix dimension mismatch");
    int m = B.dim2();

    int *indxc = new int [n];
    int *indxr = new int [n];
    int *ipiv = new int [n];

    for (int j = 0; j < n; j++)
        ipiv[j] = 0;

    for (int i = 0; i < n; i++) {
        double big = 0;
        int irow, icol;
        for (int j = 0; j < n; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < n; k++) {
                    if (ipiv[k] == 0) {
                        double Ajk = A[j][k];
                        if (std::abs(Ajk) >= big) {
                            big = std::abs(Ajk);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }
        ++ipiv[icol];

        if (irow != icol) {
            for (int j = 0; j < n; j++)
                std::swap(A[irow][j], A[icol][j]);
            for (int j = 0; j < m; j++)
                std::swap(B[irow][j], B[icol][j]);
        }

        indxr[i] = irow;
        indxc[i] = icol;
        if (A[icol][icol] == 0.0)
            error("linalg: solve_Gauss_Jordan singular matrix");
        double pivinv = 1 / A[icol][icol];
        A[icol][icol] = 1;
        for (int j = 0; j < n; j++)
            A[icol][j] *= pivinv;
        for (int j = 0; j < m; j++)
            B[icol][j] *= pivinv;

        for (int j = 0; j < n; j++) {
            if (j != icol) {
                double mult = A[j][icol];
                for (int k = 0; k < n; k++)
                    A[j][k] -= A[icol][k] * mult;
                for (int k = 0; k < m; k++)
                    B[j][k] -= B[icol][k] * mult;
            }
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        if (indxr[i] != indxc[i]) {
            for (int j = 0; j < n; j++)
                std::swap(A[j][indxr[i]], A[j][indxc[i]]);
        }
    }

    delete [] indxc;
    delete [] indxr;
    delete [] ipiv;
}

void ludcmp(Matrix<double,2>& A, int *indx, double& d) {
    const double TINY=1.0e-20;
    int i,imax,j,k;
    double big,dum,sum,temp;

    int n=A.dim1();
    Matrix<double,1> vv(n);
    d=1.0;
    for (i=0;i<n;i++) {
        big=0.0;
        for (j=0;j<n;j++)
            if ((temp=std::abs(A[i][j])) > big) big=temp;
        if (big == 0.0) error("linalg: Singular matrix in routine ludcmp");
        vv[i]=1.0/big;
    }
    for (j=0;j<n;j++) {
        for (i=0;i<j;i++) {
            sum=A[i][j];
            for (k=0;k<i;k++) sum -= A[i][k]*A[k][j];
            A[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<n;i++) {
            sum=A[i][j];
            for (k=0;k<j;k++) sum -= A[i][k]*A[k][j];
            A[i][j]=sum;
            if ((dum=vv[i]*std::abs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=0;k<n;k++) {
                dum=A[imax][k];
                A[imax][k]=A[j][k];
                A[j][k]=dum;
            }
            d = -d;
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (A[j][j] == 0.0) A[j][j]=TINY;
        if (j != n-1) {
            dum=1.0/(A[j][j]);
            for (i=j+1;i<n;i++) A[i][j] *= dum;
        }
    }
}

void lubksb(const Matrix<double,2>& A, int *indx, Matrix<double,2>& B) {
    int i,ii=0,ip,j;
    double sum;

    int n=A.dim1();
    for (int k = 0; k < B.dim2(); k++) {
        for (i=0;i<n;i++) {
            ip=indx[i];
            sum=B[ip][k];
            B[ip][k]=B[i][k];
            if (ii != 0)
                for (j=ii-1;j<i;j++) sum -= A[i][j]*B[j][k];
            else if (sum != 0.0)
                ii=i+1;
            B[i][k]=sum;
        }
        for (i=n-1;i>=0;i--) {
            sum=B[i][k];
            for (j=i+1;j<n;j++) sum -= A[i][j]*B[j][k];
            B[i][k]=sum/A[i][i];
        }
    }
}

void solve_LU_decompose(Matrix<double,2>& A, Matrix<double,2>& B) {

    int n = A.dim1();
    if (A.dim2() != n)
        error("linalg: matrix dimension mismatch");
    if (B.dim1() != n)
        error("linalg: matrix dimension mismatch");
    int m = B.dim2();

    int *indx = new int [n];
    double d;
    ludcmp(A, indx, d);
    lubksb(A, indx, B);

    delete [] indx;
}

void reduce_Householder(Matrix<double,2>& A, Matrix<double,1>& d, Matrix<double,1>& e) {

    int n = d.dim1();

    for (int i = n - 1; i > 0; i--) {
        int l = i - 1;
        double h = 0;
        double scale = 0;
        if (l > 0) {
            for (int k = 0; k <= l; k++)
                scale += std::abs(A[i][k]);
            if (scale == 0.0)
                e[i] = A[i][l];
            else {
                for (int k = 0; k <= l; k++) {
                    A[i][k] /= scale;
                    h += A[i][k] * A[i][k];
                }
                double f = A[i][l];
                double g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
                e[i] = scale * g;
                h -= f * g;
                A[i][l] = f - g;
                f = 0.0;
                for (int j = 0; j <= l; j++) {
                    A[j][i] = A[i][j] / h;
                    g = 0.0;
                    for (int k = 0; k <= j; k++)
                       g += A[j][k] * A[i][k];
                    for (int k = j + 1; k <= l; k++)
                       g += A[k][j] * A[i][k];
                    e[j] = g / h;
                    f += e[j] * A[i][j];
                }
                double hh = f / (h + h);
                for (int j = 0; j <= l; j++) {
                    f = A[i][j];
                    e[j] = g = e[j] - hh * f;
                    for (int k = 0; k <= j; k++)
                        A[j][k] -= f * e[k] + g * A[i][k];
                }
            }
        } else
            e[i] = A[i][l];
        d[i] = h;
    }
    d[0] = 0.0;
    e[0]=0.0;
    for (int i = 0; i < n; i++) {
        if (d[i] != 0.0) {
            for (int j = 0; j < i; j++) {
                double g = 0;
                for (int k = 0; k < i; k++)
                    g += A[i][k] * A[k][j];
                for (int k = 0; k < i; k++)
                    A[k][j] -= g * A[k][i];
            }
        }
        d[i] = A[i][i];
        A[i][i] = 1.0;
        for (int j = 0; j < i; j++)
            A[j][i] = A[i][j] = 0.0;
    }
}

static double pythag(double a, double b) {
    double absa = std::abs(a);
    double absb = std::abs(b);
    if (absa > absb) {
        double ratio = absb / absa;
        return absa * sqrt(1 + ratio * ratio);
    } else {
        if (absb == 0.0)
            return 0.0;
        else {
            double ratio = absa / absb;
            return absb * sqrt(1 + ratio * ratio);
        }
    }
}

inline double sign(double a, double b) {
    if (b >= 0.0) {
        if (a >= 0.0)
            return a;
        else
            return -a;
    } else {
        if (a >= 0.0)
            return -a;
        else
            return a;
    }
}

void solve_TQLI(Matrix<double,1>& d, Matrix<double,1>& e, Matrix<double,2>& z) {

    int n = d.dim1();
    for (int i = 1; i < n; i++)
        e[i-1] = e[i];
    e[n-1] = 0.0;
    for (int l = 0; l < n; l++) {
        int iter = 0, m;
        do {
            for (m = l ; m < n-1; m++) {
                double dd = std::abs(d[m]) + std::abs(d[m+1]);
                if ((std::abs(e[m]) + dd) == dd)
                    break;
            }
            if (m != l) {
                if (iter++ == 30)
                     error("linalg: Too many iterations in solve_TQLI");
                double g = (d[l+1] - d[l]) / (2.0 * e[l]);
                double r = pythag(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + sign(r, g));
                double s = 1.0;
                double c = 1.0;
                double p = 0.0;
                int i;
                for (i = m-1; i >= l; i--) {
                    double f = s * e[i];
                    double b = c * e[i];
                    e[i+1] = r = pythag(f, g);
                    if (r == 0.0) {
                        d[i+1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d[i+1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    d[i+1] = g + (p = s * r);
                    g = c * r - b;
                    for (int k = 0; k < n; k++) {
                        f = z[k][i+1];
                        z[k][i+1] = s * z[k][i] + c * f;
                        z[k][i] = c * z[k][i] - s * f;
                    }
                }
                if (r == 0.0 && i >= l)
                    continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }
}

void sort_eigen(Matrix<double,1>& d, Matrix<double,2>& z) {

    // sorts eigenvalues and eigenvector in descending order
    int n = d.dim1();
    if (z.dim1() != n || z.dim2() != n)
        error("linalg: Bad vector, matrix dimensions in sort_eigen");
    for (int i = 0; i < n - 1; i++) {
        int k = i;
        double p = d[k];
        for (int j = i; j < n; j++)
            if (d[j] >= p)
                p = d[k = j];
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (int j = 0; j < n; j++) {
                p = z[j][i];
                z[j][i] = z[j][k];
                z[j][k] = p;
            }
        }
    }
}

Matrix<double,1> solve_eigen_symmetric(Matrix<double,2>& A)
{
    // solves eigenvalue problem A x = e x
    // returns vector of eigenvalues and replaces A by matrix of eigenvectors
    int n = A.dim1();
    Matrix<double,1> d(n), e(n);
    reduce_Householder(A, d, e);
    solve_TQLI(d, e, A);
    sort_eigen(d, A);
    return d;
}

Matrix<double,1> solve_eigen_generalized(Matrix<double,2>& A, Matrix<double,2>& S)
{
    // solves generalized eigenvalue problem Ax = e S x
    int n = A.dim1();
    Matrix<double,2> V = S;
    Matrix<double,1> mu = solve_eigen_symmetric(V);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            V[i][j] /= sqrt(mu[j]);
    Matrix<double,2> VT = transpose(V);
    A = VT * A * V;
    Matrix<double,1> e = solve_eigen_symmetric(A);
    A = V * A;
    return e;
}

    inline bool is_power_of_2(int x)
    {
	    return ( (x > 0) && ((x & (x - 1)) == 0) );
    }

    static Matrix<complex<double>,1> dft(
        const Matrix<complex<double>,1>& v,
        const bool inverse=false)
    {
        int N = v.dim1();
        const double pi = 4 * atan(1.0);
        complex<double> W(0, 2 * pi / double(N));
        W = inverse ? exp(-W) : exp(W);
        Matrix<complex<double>,1> dft(N), roots_of_unity(N);
        complex<double> W_i(1, 0);
        for (int i = 0; i < N; i++) {
            roots_of_unity[i] = W_i;
            W_i *= W;
        }
        for (int k = 0; k < N; k++)
            for (int j = 0; j < N; j++)
                dft[k] += roots_of_unity[(j*k) % N] * v[j];
        return dft;
    }

    static Matrix<complex<double>,1> fft_recursive(
        const Matrix<complex<double>,1>& v,
        const bool inverse=false)
    {
        int N = v.dim1();
        if (N == 1)                 // transform is trivial
            return v;
        else if (N % 2 == 1)        // transform is odd, lemma does not apply
            return dft(v, inverse);
        // perform even-odd decomposition and transform recursively
        Matrix<complex<double>,1> even(N / 2);
        for (int j = 0; j < N / 2; j++)
            even[j] = v[2 * j];
        even = fft_recursive(even, inverse);
        Matrix<complex<double>,1> odd(N / 2);
        for (int j = 0; j < N / 2; j++)
            odd[j] = v[2 * j + 1];
        odd = fft_recursive(odd, inverse);
        const double pi = 4 * std::atan(1.0);
        const std::complex<double> i(0.0, 1.0);
        complex<double> W = 2 * pi / double(N);
        W = inverse ? exp(-i * W) : exp(i * W);
        Matrix<complex<double>,1> fft(N);
        complex<double> W_k = 1.0;
        for (int k = 0; k < N; k++) {
            fft[k] = even[k % (N / 2)] + W_k * odd[k % (N / 2)];
            W_k *= W;
        }
        if (inverse)
            fft /= 2.0;
        return fft;
    }

    static void Danielson_Lanczos(
        Matrix<complex<double>,1>& v,
        int two_pow_n,
        bool inverse)
    {
        const double pi = 4 * atan(1.0);
        int N = v.dim1();
        complex<double> W(0, pi / two_pow_n), W_j(1, 0);
        W = inverse ? exp(W) : exp(-W);
        for (int j = 0; j < two_pow_n; j++) {
            for (int i = j; i < N; i += 2 * two_pow_n) {
                complex<double> temp = W_j * v[two_pow_n + i];
                v[two_pow_n + i] = v[i] - temp;
                v[i] += temp;
            }
            W_j *= W;
        }
    }

    static void bit_reverse(
        Matrix<complex<double>,1>& v)
    {
        int N = v.dim1();
        int j = 1;
        for (int i = 1; i < N; ++i) {
            if (i < j)
                swap(v[i-1], v[j-1]);
            int k = N / 2;
            while ( k < j ) {
                j -= k;
                k /= 2;
            }
            j += k;
        }
    }

    static void fft_two_pow(
        Matrix<complex<double>,1>& v,
        bool inverse)
    {
        int N = v.dim1();
        bit_reverse(v);
        for (int two_pow_n = 1; two_pow_n < N; two_pow_n *= 2)
            Danielson_Lanczos(v, two_pow_n, inverse);
    }

    void fft(Matrix<std::complex<double>,1>& v, const bool inverse)
    {
        int N = v.dim1();
        if (N == 1)
            return;
        if (is_power_of_2(N))
            fft_two_pow(v, inverse);
        else
            v = fft_recursive(v, inverse);
        v /= sqrt(double(N));
    }

    void fft(Matrix<std::complex<double>,2>& m, const bool inverse)
    {

    }

    // Least Squares and Chi-Square Fits

    void least_squares_fit(             // makes a linear least-squares fit
        const Matrix<double,1>& x,      // vector of x values - input
        const Matrix<double,1>& y,      // vector of y values - input
        double& a,                      // fitted intercept - output
        double& b,                      // fitted slope - output
        double& sigma_a,                // estimated error in intercept - output
        double& sigma_b,                // estimated error in slope - output
        double& sigma)                  // estimated error bar in y
    {
        // determine the number of data points
        int n = x.dim1();

        // declare and initialize various sums to be computed
        double s_x = 0, s_y = 0, s_xx = 0, s_xy = 0;

        // compute the sums
        for (int i = 0; i < n; i++) {
            s_x += x[i];
            s_y += y[i];
            s_xx += x[i] * x[i];
            s_xy += x[i] * y[i];
        }

        // evaluate least-squares fit forumlas
        double denom = n * s_xx - s_x * s_x;
        a = (s_xx * s_y - s_x * s_xy) / denom;
        b = (n * s_xy - s_x * s_y) / denom;

        // estimate the variance in the data set
        double sum = 0;
        for (int i = 0; i < n; i++) {
            double y_of_x_i = a + b * x[i];
            double error = y[i] - y_of_x_i;
            sum += error * error;
        }
        sigma = sqrt(sum / (n - 2));    // estimate of error bar in y

        // estimate errors in a and b
        sigma_a = sqrt(sigma * sigma * s_xx / denom);
        sigma_b = sqrt(sigma * sigma * n / denom);
    }

    void chi_square_fit(                // makes a linear chi-square fit
        const Matrix<double,1>& x,      // vector of x values - input
        const Matrix<double,1>& y,      // vector of y values - input
        const Matrix<double,1>& err,    // vector of y error values - input
        double& a,                      // fitted intercept - output
        double& b,                      // fitted slope - output
        double& sigma_a,                // estimated error in intercept - output
        double& sigma_b,                // estimated error in slope - output
        double& chi_square)             // minimized chi-square sum - output
    {
        int n = x.dim1();

        double S = 0, S_x = 0, S_y = 0;
        for (int i = 0; i < n; i++) {
            S += 1 / err[i] / err[i];
            S_x += x[i] / err[i] / err[i];
            S_y += y[i] / err[i] / err[i];
        }

        Matrix<double,1> t(n);
        for (int i = 0; i < n; i++)
            t[i] = (x[i] - S_x/S) / err[i];

        double S_tt = 0;
        for (int i = 0; i < n; i++)
            S_tt += t[i] * t[i];

        b = 0;
        for (int i = 0; i < n; i++)
            b += t[i] * y[i] / err[i];
        b /= S_tt;

        a = (S_y - S_x * b) / S;
        sigma_a = sqrt((1 + S_x * S_x / S / S_tt) / S);
        sigma_b = sqrt(1 / S_tt);

        chi_square = 0;
        for (int i = 0; i < n; i++) {
            double diff = (y[i] - a - b * x[i]) / err[i];
            chi_square += diff * diff;
        }
    }

} /* namespace cpt */
