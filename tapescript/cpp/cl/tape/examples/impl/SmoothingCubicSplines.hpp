
#ifndef cl_tape_examples_impl_SmoothingCubicSplines_hpp
#define cl_tape_examples_impl_SmoothingCubicSplines_hpp

#define CL_BASE_SERIALIZER_OPEN
#include <cl/tape/tape.hpp>

namespace cl
{
	struct Spline{
		double a, b, c, d, x, y;
		Spline(const double& x0, const double& y0)
		{
			x = x0;
			y = y0;
		}
	};

	class CubicSplines
	{
		std::vector<Spline> S;
		std::vector<double> data_x;
		std::vector<double> data_y;
		int n;
		// Status of calculation.
		bool is_calculated = false;
	public:
		CubicSplines() = default;
		CubicSplines(const std::vector<double>& x0, const std::vector<double>& y0) : n(x0.size() - 1), data_x(x0), data_y(y0)
		{
			is_calculated = false;
			for (auto i = 0; i < n + 1; ++i) S.emplace_back(data_x[i], data_y[i]);
		}
		void load(const std::vector<double>& x0, const std::vector<double>& y0)
		{
			is_calculated = false;
			n = x0.size() - 1;
			data_x = x0;
			data_y = y0;
			for (auto i = 0; i < n + 1; ++i) S.emplace_back(data_x[i], data_y[i]);
		}
		std::vector<std::vector<double>> interpolation()
		{
			std::vector<double> h(n), p(n), q(n), b(n + 1);
			h[0] = S[1].x - S[0].x;
			for (auto i = 1; i < n; ++i)
			{
				h[i] = S[i + 1].x - S[i].x;
				p[i] = 2 * (S[i + 1].x - S[i - 1].x);
				q[i] = 3 * (S[i + 1].y - S[i].y) / h[i] - 3 * (S[i].y - S[i - 1].y) / h[i - 1];
			}
			//Gaussian Elimination
			for (auto i = 2; i < n; ++i)
			{
				p[i] -= h[i - 1] * h[i - 1] / p[i - 1];
				q[i] -= q[i - 1] * h[i - 1] / p[i - 1];
			}
			//Backsubstitution
			b[n - 1] = q[n - 1] / p[n - 1];
			for (auto i = 2; i < n; ++i) b[n - i] = (q[n - i] - h[n - i] * b[n - i + 1]) / p[n - i];
			//Spline Parameters
			S[0].a = b[1] / (3 * h[0]);
			S[0].b = 0;
			S[0].c = (S[1].y - S[0].y) / h[0] - b[1] * h[0] / 3;
			S[0].d = S[0].y;
			S[n].b = 0;
			for (auto i = 1; i < n; ++i)
			{
				S[i].a = (b[i + 1] - b[i]) / (3 * h[i]);
				S[i].b = b[i];
				S[i].c = (b[i] + b[i - 1])*h[i - 1] + S[i - 1].c;
				S[i].d = S[i].y;
			}
			std::vector<std::vector<double>> interp;
			for(auto e : S)
			{
				std::vector<double> tmp = { e.a, e.b, e.c, e.d };
				interp.push_back(tmp);
			}
			interp.pop_back();
			is_calculated = true;
			return interp;
		}
	};

	class SmoothingCubicSplines
	{
		std::vector<Spline> S;
		std::vector<double> data_x;
		std::vector<double> data_y;
		std::vector<double> sigma;
		double lambda = -1;
		int n;
		bool is_calculated = false;

		void Quincunx(std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& q)
		{
			//factorisation
			u[0] = 0;
			v[1] /= u[1];
			w[1] /= u[1];
			for (auto i = 2; i < n; ++i)
			{
				u[i] -= u[i - 2] * w[i - 2] * w[i - 2] + u[i - 1] * v[i - 1] * v[i - 1];
				v[i] = (v[i] - u[i - 1] * v[i - 1] * w[i - 1]) / u[i];
				w[i] /= u[i];
			}
			//forward substitution + maybe some actions
			for (auto i = 2; i < n; ++i) q[i] -= v[i - 1] * q[i - 1] + w[i - 2] * q[i - 2];
			for (auto i = 1; i < n; ++i) q[i] /= u[i];
			//back substitution
			q[n] = 0;
			for (auto i = n - 2; i>0; --i) q[i] -= v[i] * q[i + 1] + w[i] * q[i + 2];
		}
	public:
		
		SmoothingCubicSplines() = default;

		SmoothingCubicSplines(const std::vector<double>& x0, const std::vector<double>& y0) : n(x0.size() - 1), data_x(x0), data_y(y0)
		{
			for (auto i = 0; i < n + 1; ++i)
			{
				S.emplace_back(data_x[i], data_y[i]);
				sigma.push_back(1);
			}
			lambda = -1;
			is_calculated = false;
		}
		
		SmoothingCubicSplines(const std::vector<double>& x0, const std::vector<double>& y0, const double& lambda0) : n(x0.size() - 1), data_x(x0), data_y(y0), lambda(lambda0 >= 0 && lambda0 <= 1 ? lambda0 : -1)
		{
			for (auto i = 0; i < n + 1; ++i)
			{
				S.emplace_back(data_x[i], data_y[i]);
				sigma.push_back(1);
			}
			is_calculated = false;
		}

		void load(const std::vector<double>& x0, const std::vector<double>& y0)
		{
			reset_data();
			n = x0.size() - 1;
			data_x = x0;
			data_y = y0;
			for (auto i = 0; i < n + 1; ++i)
			{
				S.emplace_back(data_x[i], data_y[i]);
				sigma.push_back(1);
			}
			lambda = -1;
		}

		void load(const std::vector<double>& x0, const std::vector<double>& y0, const double& lambda0)
		{
			load(x0, y0);
			lambda = lambda0;
		}
		
		void smoothing_spline(const double& lambda0)
		{
			if (lambda0 >= 0 && lambda0 <= 1)
			{
				lambda = lambda0;
				smoothing_spline();
			}
		}

		void smoothing_spline()
		{
			std::vector<double> h(n), r(n + 1), f(n + 1), p(n), q(n + 1), u(n), v(n), w(n);
			double mu = lambda == 0 ? 999999999 : 2 * (1 - lambda) / (3 * lambda);
			h[0] = S[1].x - S[0].x;
			r[0] = 3 / h[0];
			for (auto i = 1; i < n; ++i)
			{
				h[i] = S[i + 1].x - S[i].x;
				r[i] = 3 / h[i];
				f[i] = -(r[i - 1] + r[i]);
				p[i] = 2 * (S[i + 1].x - S[i - 1].x);
				q[i] = 3 * (S[i + 1].y - S[i].y) / h[i] - 3 * (S[i].y - S[i - 1].y) / h[i - 1];
			}
			for (auto i = 1; i < n; ++i)
			{
				u[i] = r[i - 1] * r[i - 1] * sigma[i - 1] + f[i] * f[i] * sigma[i] + r[i] * r[i] * sigma[i + 1];
				u[i] = mu * u[i] + p[i];
				v[i] = f[i] * r[i] * sigma[i] + r[i] * f[i + 1] * sigma[i + 1];
				v[i] = mu * v[i] + h[i];
				w[i] = mu * r[i] * r[i + 1] * sigma[i + 1];
			}
			Quincunx(u, v, w, q);
			//Spline Parameters
			S[0].d = S[0].y - mu * r[0] * q[1] * sigma[0];
			S[1].d = S[1].y - mu * (f[1] * q[1] + r[1] * q[2]) * sigma[0];
			S[0].a = q[1] / (3 * h[0]);
			S[0].b = 0;
			S[0].c = (S[1].d - S[0].d) / h[0] - q[1] * h[0] / 3;
			r[0] = 0;
			for (auto i = 1; i < n; ++i)
			{
				S[i].a = (q[i + 1] - q[i]) / (3 * h[i]);
				S[i].b = q[i];
				S[i].c = (q[i] + q[i - 1]) * h[i - 1] + S[i - 1].c;
				S[i].d = r[i - 1] * q[i - 1] + f[i] * q[i] + r[i] * q[i + 1];
				S[i].d = S[i].y - mu * S[i].d * sigma[i];
			}
			is_calculated = true;
		}

		std::vector<double> get_fitted_value()
		{
			if (is_calculated)
			{
				std::vector<double> smooth;
				for (auto e : S) smooth.push_back(e.d);
				smooth.pop_back();
				smooth.push_back(S[n - 1].a * (data_x[n] - data_x[n - 1]) * (data_x[n] - data_x[n - 1]) * (data_x[n] - data_x[n - 1]) + S[n - 1].b * (data_x[n] - data_x[n - 1]) * (data_x[n] - data_x[n - 1]) + S[n - 1].c * (data_x[n] - data_x[n - 1]) + S[n - 1 ].d);
				return smooth;
			}
			else
			{
				std::vector<double> err;
				return err;
			};
		}

		void reset_data()
		{
			S.clear();
			data_x.clear();
			data_y.clear();
			sigma.clear();
			lambda = -1;
			n = 0;
			is_calculated = false;
		}
	};
}

#endif // cl_tape_examples_impl_SmoothingCubicSplines_hpp
