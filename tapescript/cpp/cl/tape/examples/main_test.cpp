#include <iostream>

#define CL_BASE_SERIALIZER_OPEN
#define CL_USE_TAPE_SERIALIZATION

#include <cl/tape/impl/detail/enable_ad.hpp>
#include <cl/tape/tape.hpp>
#include "impl/utils.hpp"

#pragma region myFuncs
// define y(x) = Poly(a, x) in the empty namespace
template <class T>
T Poly(const std::vector<double> &a, const T &x)
{
	size_t k = a.size();
	T y = 0.;  // initialize summation
	T x_i = 1.;  // initialize x^i
	size_t i;
	for (i = 0; i < k; i++)
	{
		y += a[i] * x_i;  // y   = y + a_i * x^i
		x_i *= x;           // x_i = x_i * x
	}
	return y;
}

template <class Type>
Type exp_2(const Type &x)
{
	Type v1 = x;                // v1 = x
	Type v2 = Type(1) + v1;     // v2 = 1 + x
	Type v3 = v1 * v1;          // v3 = x^2
	Type v4 = v3 / Type(2);     // v4 = x^2 / 2
	Type v5 = v2 + v4;          // v5 = 1 + x + x^2 / 2
	return v5;                   // exp_2(x) = 1 + x + x^2 / 2
}

template <class T>
std::vector<T> MyFunc(const std::vector<T> &x)
{
	std::vector<T> y(x.size(), 0.);  // initialize summation
	for (auto i = 0; i < x.size(); ++i)
	{
		y[i] += 4 * x[i] + 7 * x[i] * x[i] * x[i]; // exp(x) - 1 - sqrt(x / (1 + x)) - x*x / T(2);
	}
	return y;
}

bool exp_2_cppad(void)
{
	bool ok = true;
	using CppAD::AD;
	using CppAD::NearEqual; // checks if values are nearly equal

	// domain space vector
	size_t n = 1; // dimension of the domain space
	std::vector< AD<double> > X(n);
	X[0] = .5;    // value of x for this operation sequence

	// declare independent variables and start recording operation sequence
	CppAD::Independent(X);

	// evaluate our exponential approximation
	AD<double> x = X[0];
	AD<double> apx = exp_2(x);

	// range space vector
	size_t m = 1;  // dimension of the range space
	std::vector< AD<double> > Y(m);
	Y[0] = apx;    // variable that represents only range space component

	// Create f: X -> Y corresponding to this operation sequence
	// and stop recording. This also executes a zero order forward
	// sweep using values in X for x.
	CppAD::ADFun<double> f(X, Y);

	// first order forward sweep that computes
	// partial of exp_2(x) with respect to x
	std::vector<double> dx(n);  // differential in domain space
	std::vector<double> dy(m);  // differential in range space
	dx[0] = 1.;            // direction for partial derivative
	dy = f.Forward(1, dx);
	double check = 1.5;
	ok &= NearEqual(dy[0], check, 1e-10, 1e-10);

	// first order reverse sweep that computes the derivative
	std::vector<double>  w(m);   // weights for components of the range
	std::vector<double> dw(n);   // derivative of the weighted function
	w[0] = 1.;              // there is only one weight
	dw = f.Reverse(1, w); // derivative of w[0] * exp_2(x)
	check = 1.5;            // partial of exp_2(x) with respect to x
	ok &= NearEqual(dw[0], check, 1e-10, 1e-10);

	// second order forward sweep that computes
	// second partial of exp_2(x) with respect to x
	std::vector<double> x2(n);     // second order Taylor coefficients
	std::vector<double> y2(m);
	x2[0] = 0.;               // evaluate second partial .w.r.t. x
	y2 = f.Forward(2, x2);
	check = 0.5 * 1.;         // Taylor coef is 1/2 second derivative
	ok &= NearEqual(y2[0], check, 1e-10, 1e-10);

	// second order reverse sweep that computes
	// derivative of partial of exp_2(x) w.r.t. x
	dw.resize(2 * n);         // space for first and second derivatives
	dw = f.Reverse(2, w);
	check = 1.;               // result should be second derivative
	ok &= NearEqual(dw[0 * 2 + 1], check, 1e-10, 1e-10);

	return ok;
}

#pragma endregion

int main()
{
	auto cl_start = clock();
#pragma region get_started
	//{
	//	using CppAD::AD;           // use AD as abbreviation for CppAD::AD
	//	using std::vector;         // use vector as abbreviation for vector
	//	size_t i;                  // a temporary index

	//	// vector of polynomial coefficients
	//	size_t k = 5;              // number of polynomial coefficients
	//	vector<double> a(k);       // vector of polynomial coefficients
	//	for (i = 0; i < k; i++)
	//		a[i] = 1.;           // value of polynomial coefficients

	//	// domain space vector
	//	size_t n = 1;              // number of domain space variables
	//	vector< AD<double> > X(n); // vector of domain space variables
	//	X[0] = 3.;                 // value corresponding to operation sequence

	//	// declare independent variables and start recording operation sequence
	//	CppAD::Independent(X);

	//	// range space vector
	//	size_t m = 1;              // number of ranges space variables
	//	vector< AD<double> > Y(m); // vector of ranges space variables
	//	Y[0] = Poly(a, X[0]);      // value during recording of operations

	//	// store operation sequence in f: X -> Y and stop recording
	//	CppAD::ADFun<double> f(X, Y);

	//	// compute derivative using operation sequence stored in f
	//	vector<double> jac(m * n); // Jacobian of f (m by n matrix)
	//	vector<double> x(n);       // domain space vector
	//	x[0] = 3.;                 // argument value for derivative
	//	jac = f.Jacobian(x);      // Jacobian for operation sequence

	//	// print the results
	//	std::cout << "f'(" << x[0] << ") computed by CppAD = " << jac[0] << std::endl;

	//	// check if the derivative is correct
	//	int error_code;
	//	if (jac[0] == 142.)
	//		error_code = 0;      // return code for correct case
	//	else  error_code = 1;      // return code for incorrect case
	//}
#pragma endregion

#pragma region get_started_my
	{
		using CppAD::AD;           // use AD as abbreviation for CppAD::AD

		// domain space vector
		size_t n = 1;              // number of domain space variables
		std::vector< AD<double> > X(n); // vector of domain space variables
		X[0] = 5.;                 // value corresponding to operation sequence

		// declare independent variables and start recording operation sequence
		CppAD::Independent(X);

		// range space vector
		size_t m = 1;              // number of ranges space variables
		std::vector< AD<double> > Y(m); // vector of ranges space variables
		Y = MyFunc(X); // value during recording of operations

		// store operation sequence in f: X -> Y and stop recording
		CppAD::ADFun<double> f(X, Y);

		// compute derivative using operation sequence stored in f
		std::vector<double> jac(m * n); // Jacobian of f (m by n matrix)
		std::vector<double> x(n);       // domain space vector
		x[0] = 1.;                 // argument value for derivative
		jac = f.Forward(1, x); // .Jacobian(x);      // Jacobian for operation sequence

		// print the results
		std::cout << "MyFunc'(" << X[0] << ") computed by CppAD = " << jac[0] << std::endl;
	}
#pragma endregion

#pragma region test_exp_2
	//{
	//	bool ok = true;

	//	ok &= Run(exp_2_cppad, "exp_2_cppad");
	//	if (ok) cout << "All " << int(Run_ok_count) << " tests passed." << endl;
	//	else cout << int(Run_error_count) << " tests failed." << endl;
	//}
#pragma endregion

#pragma region tmp
	{

	}
#pragma endregion

	std::cout << "\nTime: " << clock() - cl_start << " ms" << std::endl;
	system("pause");
	return 0;
}
