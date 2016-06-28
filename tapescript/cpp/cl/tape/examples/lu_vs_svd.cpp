#include <iostream>

#define CL_BASE_SERIALIZER_OPEN
#define CL_USE_TAPE_SERIALIZATION

#include <cl/tape/impl/detail/enable_ad.hpp>
#include <cl/tape/tape.hpp>
#include "impl/utils.hpp"
#include "impl/polynomial_regression_examples.hpp"

namespace cl
{
	// Invert symmetric matrix using LU decomposition from the Boost library.
	// This template method works with AD and not AD value_type. // NO, just AD coz for not AD are 2 overloads operator '<'
	template<typename Matrix>
	static Matrix invert_sym_matrix(const Matrix& mat)
	{
		typedef typename Matrix::value_type Array;
		typedef typename Array::value_type element_type;
		// Matrix size.
		int m = mat.size();
		boost::numeric::ublas::matrix<element_type> input(m, m);
		for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
			input(i, j) = mat[i][j];
		// Create a permutation matrix for the LU-factorization
		boost::numeric::ublas::permutation_matrix<std::size_t> pm(m);
		// Perform LU-factorization
		if (boost::numeric::ublas::lu_factorize(input, pm))
			throw std::runtime_error("Singular matrix");
		// Create identity matrix of for inverse
		boost::numeric::ublas::matrix<element_type> inverse(m, m);
		inverse.assign(boost::numeric::ublas::identity_matrix<element_type>(m));
		try
		{
			boost::numeric::ublas::lu_substitute(input, pm, inverse);
		}
		catch (boost::numeric::ublas::internal_logic&)
		{
			// Singular matrix may not be detected in lu_factorize due to the round error.
			// In that case lu_substitute throws internal_logic exception.
			throw "internal_logic error";
		}
		Matrix mat_inv = mat;
		for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
			mat_inv[i][j] = inverse(i, j);
		return mat_inv;
	}

	// Invert symmetric matrix using SV decomposition.
	// This template method works with AD and not AD value_type.
	template<typename Matrix>
	static Matrix invert_sym_matrix_SVD(const Matrix& mat)
	{
		typedef typename Matrix::value_type Array;
		typedef typename Array::value_type element_type;
		SVD<element_type> decomp(mat);
		return decomp.invert_matrix();
	}
}

// Reads data from .csv file.
template<typename Matrix>
Matrix read(std::string name)
{
	typedef typename Matrix::value_type Array;
	typedef typename Array::value_type T;
	Matrix values;
	Array valueline;
	std::ifstream fin("input/" + name + ".csv");
	std::string item;
	for (std::string line; getline(fin, line);)
	{
		std::istringstream in(line);
		while (getline(in, item, ','))
			valueline.push_back(atof(item.c_str()));
		values.push_back(valueline);
		valueline.clear();
	}
	return values;
}

// Formatted matrix output.
template<class T>
void printMatrix(const boost::numeric::ublas::matrix<T> &m, std::ostream & out = std::cout)
{
	for (unsigned i = 0; i < m.size1(); ++i)
	{
		out << "| ";
		for (unsigned j = 0; j < m.size2(); ++j)
		{
			out << std::setw(14) << std::setprecision(7) << m(i, j) << " | ";
		}
		out << "|" << std::endl;
	}
	out << std::endl;
}

#pragma region SVD
template<class T>
class SVD
{
public:
	int row_number, col_number;
	// matrices of left singular std::vectors and right singular std::vectors.
	std::vector<std::vector<T>> u, v;
	// std::vector of singular values.
	std::vector<T> w;
	// epsilon and treshhold.
	T eps, tsh;
	// constructor with initial matrix a. Initializes matrices, makes the decompose and the reorder.
	SVD(const std::vector<std::vector<T>> &a) : row_number(a.size()), col_number(a[0].size()), u(a), v(col_number), w(col_number)
	{
		for (int i = 0; i < col_number; i++)
			v[i].resize(col_number);
		eps = std::numeric_limits<T>::epsilon();
		decompose();
		reorder();
		// the treshhold = 1 / 2 * eps * w(1) * sqrt(n + m + 1) where w(1) is the biggest singular value; n, m - the number of rows and columns in the matrix a; eps - very small value.
		tsh = 0.5 * sqrt(row_number + col_number + 1.) * w[0] * eps;
	}
	// solves the linear system A * x = b, where x and b are std::vectors.
	void solve(std::vector<T> &b, std::vector<T> &x, T thresh);
	// solves the linear system A * x = b, where x and b are matrices.
	void solve(std::vector<std::vector<T>> &b, std::vector<std::vector<T>> &x, T thresh);
	// invert A to A-1
	std::vector<std::vector<T>> invert_matrix();
	// returns the rank of the matrix.
	int rank(T thresh);
	// returns the number of singular values that less than treshhold.
	int nullity(T thresh);
	// returns the matrix where singular values are more than the treshhold.
	std::vector<std::vector<T>> range(T thresh);
	// returns the matrix where singular values are more than the treshhold.
	std::vector<std::vector<T>> nullspace(T thresh);
	// returns the ration between the first and the last singular value.
	T inv_condition() { return (w[0] <= 0. || w[col_number - 1] <= 0.) ? 0. : w[col_number - 1] / w[0]; }
	// makes the singular value decomposition of the matrix a.
	void decompose();
	// reorders in descent way the singular values.
	void reorder();
	// Computes sqrt(a^2 + b^2) without underflow or overflow.
	T pythag(const T a, const T b);
};

template<class T>
inline T SIGN(const T &a, const T &b);

// solves the linear system A * x = b, where x and b are std::vectors.
template<class T>
void SVD<T>::solve(std::vector<T> &b, std::vector<T> &x, T thresh = -1.)
{
	int i, j, jj;
	T s;
	if (b.size() != row_number || x.size() != col_number) throw("SVD::solve bad sizes");
	std::vector<T> tmp(col_number);
	tsh = (thresh >= 0. ? thresh : 0.5 * sqrt(row_number + col_number + 1.) * w[0] * eps);
	for (j = 0; j < col_number; j++)
	{
		s = 0.0;
		if (w[j] > tsh)
		{
			for (i = 0; i < row_number; i++)
				s += u[i][j] * b[i];
			s /= w[j];
		}
		tmp[j] = s;
	}
	for (j = 0; j < col_number; j++)
	{
		s = 0.0;
		for (jj = 0; jj < col_number; jj++)
			s += v[j][jj] * tmp[jj];
		x[j] = s;
	}
}

// solves the linear system A * x = b, where x and b are matrices.
template<class T>
void SVD<T>::solve(std::vector<std::vector<T>> &b, std::vector<std::vector<T>> &x, T thresh = -1.)
{
	int i, j, p = b[0].size();
	if (b.size() != row_number || x[0].size() != col_number || x[0].size() != p)
		throw("SVD::solve bad sizes");
	std::vector<T> xx(col_number), bcol(row_number);
	for (j = 0; j < p; j++)
	{
		for (i = 0; i < row_number; i++)
			bcol[i] = b[i][j];
		solve(bcol, xx, thresh);
		for (i = 0; i < col_number; i++)
			x[i][j] = xx[i];
	}
}

// invert A to A-1
template<class T>
std::vector<std::vector<T>> SVD<T>::invert_matrix()
{
	std::vector<std::vector<T>> res(u);
	for (int j = 0; j < col_number; j++)
	{
		T s;
		std::vector<T> tmp(col_number);
		for (int jj = 0; jj < col_number; jj++)
		{
			s = 0.0;
			for (int ii = 0; ii < row_number; ii++)
				s += (ii == j ? u[ii][jj] : 0); // only if on diag.
			s /= w[jj];
			tmp[jj] = s;
		}
		for (int jj = 0; jj < col_number; jj++)
		{
			s = 0.0;
			for (int ii = 0; ii < col_number; ii++)
				s += v[jj][ii] * tmp[ii];
			res[jj][j] = s;
		}
	}
	return res;
}

// returns the rank of the matrix.
template<class T>
int SVD<T>::rank(T thresh = -1.)
{
	int j, nr = 0;
	tsh = (thresh >= 0. ? thresh : 0.5 * sqrt(row_number + col_number + 1.) * w[0] * eps);
	for (j = 0; j < col_number; j++)
	if (w[j] > tsh)
		nr++;
	return nr;
}

// returns the number of singular values that less than treshhold
template<class T>
int SVD<T>::nullity(T thresh = -1.)
{
	int j, nn = 0;
	tsh = (thresh >= 0. ? thresh : 0.5 * sqrt(row_number + col_number + 1.) * w[0] * eps);
	for (j = 0; j < col_number; j++)
	if (w[j] <= tsh)
		nn++;
	return nn;
}

// returns the matrix where singular values are more than the treshhold.
template<class T>
std::vector<std::vector<T>> SVD<T>::range(T thresh = -1.)
{
	int i, j, nr = 0;
	std::vector<std::vector<T>> rnge(row_number);
	for (int i = 0; i < row_number; i++)
		rnge[i].resize(rank(thresh));
	for (j = 0; j < col_number; j++)
	{
		if (w[j] > tsh)
		{
			for (i = 0; i<row_number; i++)
				rnge[i][nr] = u[i][j];
			nr++;
		}
	}
	return rnge;
}

// returns the matrix where singular values are less than the treshhold.
template<class T>
std::vector<std::vector<T>> SVD<T>::nullspace(T thresh = -1.)
{
	int j, jj, nn = 0;
	std::vector<std::vector<T>> nullsp(col_number);
	for (int i = 0; i < col_number; i++)
		nullsp.resize(nullity(thresh));
	for (j = 0; j < col_number; j++)
	{
		if (w[j] <= tsh)
		{
			for (jj = 0; jj < col_number; jj++)
				nullsp[jj][nn] = v[jj][j];
			nn++;
		}
	}
	return nullsp;
}

template<class T>
void SVD<T>::decompose()
{
	bool flag;
	int i, its, j, jj, k, l, nm;
	T anorm, c, f, g, h, s, scale, x, y, z;
	std::vector<T> rv1(col_number);
	g = scale = anorm = 0.0;
	for (i = 0; i < col_number; i++)                                                                // Householder reduction to bidiagonal form.
	{
		l = i + 2;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if (i < row_number)
		{
			for (k = i; k < row_number; k++)
				scale += abs(static_cast<double>(u[k][i]));
			if (scale != 0.0)
			{
				for (k = i; k < row_number; k++)
				{
					u[k][i] /= scale;
					s += u[k][i] * u[k][i];
				}
				f = u[i][i];
				g = -SIGN(static_cast<T>(sqrt(static_cast<double>(s))), f);
				h = f * g - s;
				u[i][i] = f - g;
				for (j = l - 1; j < col_number; j++)
				{
					for (s = 0.0, k = i; k < row_number; k++)
						s += u[k][i] * u[k][j];
					f = s / h;
					for (k = i; k < row_number; k++)
						u[k][j] += f * u[k][i];
				}
				for (k = i; k < row_number; k++)
					u[k][i] *= scale;
			}
		}
		w[i] = scale *g;
		g = s = scale = 0.0;
		if (i + 1 <= row_number && i + 1 != col_number)
		{
			for (k = l - 1; k < col_number; k++)
				scale += abs(static_cast<double>(u[i][k]));
			if (scale != 0.0)
			{
				for (k = l - 1; k < col_number; k++)
				{
					u[i][k] /= scale;
					s += u[i][k] * u[i][k];
				}
				f = u[i][l - 1];
				g = -SIGN(static_cast<T>(sqrt(static_cast<double>(s))), f);
				h = f * g - s;
				u[i][l - 1] = f - g;
				for (k = l - 1; k < col_number; k++)
					rv1[k] = u[i][k] / h;
				for (j = l - 1; j < row_number; j++)
				{
					for (s = 0.0, k = l - 1; k < col_number; k++)
						s += u[j][k] * u[i][k];
					for (k = l - 1; k < col_number; k++)
						u[j][k] += s*rv1[k];
				}
				for (k = l - 1; k < col_number; k++)
					u[i][k] *= scale;
			}
		}
		anorm = std::max(anorm, (abs(static_cast<double>(w[i])) + abs(static_cast<double>(rv1[i]))));
	}
	for (i = col_number - 1; i >= 0; i--)                                                    // Accumulation of right-hand transformations.
	{
		if (i < col_number - 1)
		{
			if (g != 0.0)
			{
				for (j = l; j < col_number; j++)
					v[j][i] = (u[i][j] / u[i][l]) / g;
				for (j = l; j < col_number; j++)
				{
					for (s = 0.0, k = l; k < col_number; k++)
						s += u[i][k] * v[k][j];
					for (k = l; k < col_number; k++)
						v[k][j] += s * v[k][i];
				}
			}
			for (j = l; j < col_number; j++)
				v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for (i = std::min(row_number, col_number) - 1; i >= 0; i--)                // Accumulation of left-hand transformations.
	{
		l = i + 1;
		g = w[i];
		for (j = l; j < col_number; j++)
			u[i][j] = 0.0;
		if (g != 0.0)
		{
			g = 1.0 / g;
			for (j = l; j < col_number; j++)
			{
				for (s = 0.0, k = l; k < row_number; k++)
					s += u[k][i] * u[k][j];
				f = (s / u[i][i])*g;
				for (k = i; k < row_number; k++)
					u[k][j] += f * u[k][i];
			}
			for (j = i; j < row_number; j++)
				u[j][i] *= g;
		}
		else for (j = i; j < row_number; j++)
			u[j][i] = 0.0;
		u[i][i] += 1;
	}
	for (k = col_number - 1; k >= 0; k--)                                                    // Diagonalization of the bidiagonal form.
	{
		for (its = 0; its < 30; its++)
		{
			flag = true;
			for (l = k; l >= 0; l--)
			{
				nm = l - 1;
				if (l == 0 || abs(static_cast<double>(rv1[l])) <= eps * anorm)                // Test for splitting.
				{
					flag = false;
					break;
				}
				if (abs(static_cast<double>(w[nm])) <= eps * anorm) break;
			}
			if (flag)                                                                                       // Cancellation of rv1[l], if l > 0.
			{
				c = 0.0;
				s = 1.0;
				for (i = l; i < k + 1; i++)
				{
					f = s * rv1[i];
					rv1[i] = c * rv1[i];
					if (abs(static_cast<double>(f)) <= eps * anorm)
						break;
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (j = 0; j < row_number; j++)
					{
						y = u[j][nm];
						z = u[j][i];
						u[j][nm] = y * c + z * s;
						u[j][i] = z * c - y * s;
					}
				}
			}
			z = w[k];
			if (l == k)                                                                                     // Convergence.
			{
				if (z < 0.0)
				{
					w[k] = -z;                                                                        // Singular value is made nonnegative.
					for (j = 0; j < col_number; j++)
						v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 29) throw("no convergence in 30 svdcmp iterations");
			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0;                                                                              // Next QR transformation
			for (j = l; j <= nm; j++)
			{
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for (jj = 0; jj < col_number; jj++)
				{
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x * c + z * s;
					v[jj][i] = z * c - x * s;
				}
				z = pythag(f, h);
				w[j] = z;
				if (z != 0)                                                                              // Rotation can be arbitrary if z D 0.
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for (jj = 0; jj < row_number; jj++)
				{
					y = u[jj][j];
					z = u[jj][i];
					u[jj][j] = y * c + z * s;
					u[jj][i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
}

template<class T>
void SVD<T>::reorder()
{
	int i, j, k, s, inc = 1;
	T sw;
	std::vector<T> su(row_number), sv(col_number);
	do
	{
		inc *= 3;
		inc++;
	}
	while (inc <= col_number);
	do
	{
		inc /= 3;
		for (i = inc; i < col_number; i++)
		{
			sw = w[i];
			for (k = 0; k < row_number; k++)
				su[k] = u[k][i];
			for (k = 0; k < col_number; k++)
				sv[k] = v[k][i];
			j = i;
			while (w[j - inc] < sw)
			{
				w[j] = w[j - inc];
				for (k = 0; k < row_number; k++)
					u[k][j] = u[k][j - inc];
				for (k = 0; k < col_number; k++)
					v[k][j] = v[k][j - inc];
				j -= inc;
				if (j < inc)
					break;
			}
			w[j] = sw;
			for (k = 0; k < row_number; k++)
				u[k][j] = su[k];
			for (k = 0; k < col_number; k++)
				v[k][j] = sv[k];

		}
	}
	while (inc > 1);
	for (k = 0; k < col_number; k++)
	{
		s = 0;
		for (i = 0; i < row_number; i++)
		if (u[i][k] < 0.)
			s++;
		for (j = 0; j < col_number; j++)
		if (v[j][k] < 0.)
			s++;
		if (s > (row_number + col_number) / 2)
		{
			for (i = 0; i < row_number; i++)
				u[i][k] = -u[i][k];
			for (j = 0; j < col_number; j++)
				v[j][k] = -v[j][k];
		}
	}
}

template<class T>
T SVD<T>::pythag(const T a, const T b)
{
	T absa = abs(static_cast<double>(a)), absb = abs(static_cast<double>(b));
	return (absa > absb ? absa * sqrt(1.0 + pow(static_cast<double>(absb / absa), 2)) :
			(absb == 0.0 ? 0.0 : absb * sqrt(1.0 + pow(static_cast<double>(absa / absb), 2))));
}

template<class T>
inline T SIGN(const T &a, const T &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
#pragma endregion

template<class T>
void SVD_test()
{
	namespace ublas = boost::numeric::ublas;
	bool release = false; // check mode for not working LU in Debug
	std::vector<std::vector<T>> a = read<std::vector<std::vector<T>>>("matrix"); // our symetric matrix
	std::vector<std::vector<T>> a1SVD, a1LU; // A-1

	// A-1 calc
	a1SVD = cl::invert_sym_matrix_SVD(a);
	if (release) a1LU = cl::invert_sym_matrix(a); // comment if T is not AD type

	// Transform to Boost matrix
	int n = a.size();
	ublas::matrix<T> a1SVDb(n, n), a1LUb(n, n), ab(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			a1SVDb(i, j) = a1SVD[i][j];
			if (release) a1LUb(i, j) = a1LU[i][j];
			ab(i, j) = a[i][j];
		}

	// A * A-1
	ublas::matrix<T> ESVDb = ublas::prod(ab, a1SVDb), ELUb;
	if (release) ELUb = ublas::prod(ab, a1LUb);

	std::ofstream out("output/lu_vs_svd1.txt");

	out << "Start A matrix:" << std::endl;
	printMatrix(ab, out);
	out << "A-1 matrix calc by SVD:" << std::endl;
	printMatrix(a1SVDb, out);
	out << "A-1 matrix calc by LU decomposition:" << std::endl;
	if (release) printMatrix(a1LUb, out);
	out << "A * A-1 matrix calc by SVD:" << std::endl;
	printMatrix(ESVDb, out);
	out << "A * A-1 matrix calc by LU decomposition:" << std::endl;
	if (release) printMatrix(ELUb, out);

	// Speed test
	auto cl_start = clock();
	long count = 1000;
	for (int i = 0; i < count; ++i)
		a1SVD = cl::invert_sym_matrix_SVD(a);
	out << "Time to calc A-1 by SVD in " << (release ? "Release" : "Debug") << " mode (" << count << " times) : " << clock() - cl_start << " ticks" << std::endl;
	cl_start = clock();
	if (release)
	{
		for (int i = 0; i < count; ++i)
			a1LU = cl::invert_sym_matrix(a); // comment if T is not AD type
		out << "Time to calc A-1 by LU decomposition in Release mode (" << count << " times) : " << clock() - cl_start << " ticks" << std::endl;
	}
	out.close();
}

int main()
{
	SVD_test<cl::tdouble>();

	system("pause");
	return 0;
}
