/*
Copyright (C) 2003-2015 CompatibL

This file is part of TapeScript, an open source library and tape encoding
standard for adjoint algorithmic differentiation (AAD), available from

http://github.com/compatibl/tapescript (source)
http://tapescript.org (documentation)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef cl_tape_examples_impl_cubic_spline_utils_hpp
#define cl_tape_examples_impl_cubic_spline_utils_hpp

#include <cl/tape/tape.hpp>
#include <cl/tape/impl/math/cubic_spline_math.hpp>
#include "impl/utils.hpp"

namespace cl
{
    class smooth_cubic_spline_utils
    {
    protected:
        // Number of input points.
        int data_number;
        // Spline sections.
        std::vector<tobject> x;
		std::vector<tobject> y;
		// u, v and w params
		std::vector<std::vector<tobject>> uvw_;
        // Auxiliary variables.
        std::vector<std::vector<tobject>> auxiliaries;
        // Spline coefficients.
        std::vector<std::vector<tobject>> coefficients;
        // Fitted values.
        std::vector<tobject> y_fitted;

    public:

        smooth_cubic_spline_utils(){}
        virtual ~smooth_cubic_spline_utils(){}

        virtual std::vector<tobject> get_fitted_values()
        {
            return y_fitted;
        }

        // Spline initialisation from X, Y input arrays.
        virtual void load(const std::vector<tobject>& x0, const std::vector<tobject>& y0)
        {
            // Number of points.
            data_number = x0.size();

            // Initialise spline input data.
            x = x0;
            y = y0;
        }

        // Calculate spline coefficients.
        virtual void Solve(tobject lambda)
        {
            coefficients.resize(3);
            for (int i = 0; i < 3; i++)
				coefficients[i].resize(data_number - 1); // n - 1 -> n
            // Calculate spline helper variables.
			calculate_helper_variables(x, auxiliaries);
            calculate_q();
            // Bound smoothing parameter value to avoid problems with NaN.
            tobject lambdaMin = 1e-7;
            tobject lambdaMax = 1 - 1e-7;
            if (lambda < lambdaMin)
                lambda = lambdaMin;
            else if (lambda > lambdaMax)
                lambda = lambdaMax;

            // Internal smoothing parameter: transform lambda -> mu = 2 * (1 - lambda) / (3 * lambda).
            tobject mu = 2 * (1 - lambda) / (3 * lambda);

            // Prepare system of linear equations: calculate diagonal, first and second upper diagonal vectors u, v, w, respectively.
            prepare_linear_system(mu);

            // Find b-coefficients in splines in auxiliaries[4] array.
            cl::tapescript::find_b_coefficients<tobject>(uvw_, auxiliaries);
            std::vector<tobject> d_coeffs;
            // Find the rest of spline parameters.
            find_all_coefficients(mu, d_coeffs);
            coefficients.push_back(d_coeffs);
            fitted_values();
        }

        virtual void fitted_values()
        {
            int n = data_number;
            y_fitted.resize(n);
            for (int i = 0; i < n; i++)
                y_fitted[i] = value_at(x[i]);
        }

        /// <summary>Get spline derivative at point x with possible extrapolation.</summary>
		virtual tobject value_at(tobject point)
        {
            int n = data_number;
            if (point < x[0])
            {
                tobject diff = point - x[0];
                return coefficients[0][0] * std::pow(diff, 3) + coefficients[1][0] * std::pow(diff, 2) +
                    coefficients[2][0] * diff + coefficients[3][0];
            }
            else if (point > x[n - 1])
            {
                tobject diff = point - x[n - 1];
                return coefficients[0][n - 2] * std::pow(diff, 3) + coefficients[1][n - 2] * std::pow(diff, 2) +
                    coefficients[2][n - 2] * diff + coefficients[3][n - 2];
            }

            // Check for requested x within spline.
            for (int j = 0; j < data_number - 1; j++)
                if (point >= x[j] && point <= x[j + 1])
                    return value_in_section(point, j);

            // If spline section for requested x was not found, return empty.
            throw std::exception("Value not found.");
        }

		virtual tobject value_in_section(tobject point, int section)
        {
            tobject diff = point - x[section];
            return coefficients[0][section] * std::pow(diff, 3) + coefficients[1][section] * std::pow(diff, 2) +
                coefficients[2][section] * diff + coefficients[3][section];
        }

		//!!!!!!!!!!!!!!!!!!!!!!!!//
		virtual void calculate_q()
		{
			int n = y.size();
			for (int i = 1; i < n - 1; ++i)
				auxiliaries[4][i - 1] = 3 * (y[i + 1] - y[i]) / auxiliaries[1][i] - 3 * (y[i] - y[i - 1]) / auxiliaries[1][i - 1];
		}
		// Prepare system of linear equations: calculate all non-zero matrix elements.
		virtual void prepare_linear_system(tobject mu)
		{
			int n = auxiliaries[0].size() + 2;
			uvw_.resize(3);
			uvw_[0].resize(n - 2);
			uvw_[1].resize(n - 3);
			uvw_[2].resize(n - 4);
			std::vector<tobject>& u = uvw_[0];
			std::vector<tobject>& v = uvw_[1];
			std::vector<tobject>& w = uvw_[2];
			std::vector<tobject>& p = auxiliaries[0];
			std::vector<tobject>& h = auxiliaries[1];
			std::vector<tobject>& r = auxiliaries[2];
			std::vector<tobject>& f = auxiliaries[3];
			for (int i = 0; i < n - 4; i++)
			{
				u[i] = mu * (std::pow(r[i], 2) + std::pow(f[i], 2) + std::pow(r[i + 1], 2)) + p[i];
				v[i] = mu * r[i + 1] * (f[i] + f[i + 1]) + h[i + 1];
				w[i] = mu * r[i + 1] * r[i + 2];
			}
			u[n - 4] = mu * (std::pow(r[n - 4], 2) + std::pow(f[n - 4], 2) + std::pow(r[n - 3], 2)) + p[n - 4];
			u[n - 3] = mu * (std::pow(r[n - 3], 2) + std::pow(f[n - 3], 2) + std::pow(r[n - 2], 2)) + p[n - 3];
			v[n - 4] = mu * r[n - 3] * (f[n - 4] + f[n - 3]) + h[n - 3];
		}
		// Calculate all spline coefficients.
		virtual void find_all_coefficients(tobject mu, std::vector<tobject>& D)
		{
			D = cl::tapescript::find_d_coefficients<std::vector<tobject>, tobject>(mu, y, auxiliaries);
			// Compact variable name for number of data points.
			int n = D.size() + 1;
			coefficients.resize(3);
			for (int i = 0; i < 3; i++)
				coefficients[i].resize(n - 1);
			coefficients[0][0] = auxiliaries[4][0] / (3 * auxiliaries[1][0]);
			coefficients[1][0] = 0.0;
			coefficients[2][0] = (D[1] - D[0]) / auxiliaries[1][0] - auxiliaries[4][0] * auxiliaries[1][0] / 3;
			coefficients[0][1] = (auxiliaries[4][1] - auxiliaries[4][0]) / (3 * auxiliaries[1][1]);
			coefficients[1][1] = auxiliaries[4][0];
			coefficients[2][1] = auxiliaries[4][0] * auxiliaries[1][0] + coefficients[2][0];
			for (int i = 2; i < n - 2; i++)
			{
				coefficients[0][i] = (auxiliaries[4][i] - auxiliaries[4][i - 1]) / (3 * auxiliaries[1][i]);
				coefficients[1][i] = auxiliaries[4][i - 1];
				coefficients[2][i] = (auxiliaries[4][i - 1] + auxiliaries[4][i - 2]) * auxiliaries[1][i - 1] + coefficients[2][i - 1];
			}
			// Since q[n - 1] = 0.
			coefficients[0][n - 2] = (0 - auxiliaries[4][n - 3]) / (3 * auxiliaries[1][n - 2]);
			coefficients[1][n - 2] = auxiliaries[4][n - 3];
			coefficients[2][n - 2] = (auxiliaries[4][n - 3] + auxiliaries[4][n - 4]) * auxiliaries[1][n - 3] + coefficients[2][n - 3];
		}

		virtual void calculate_helper_variables(const std::vector<tobject>& x0, std::vector<std::vector<tobject>>& aux)
		{
			cl::tapescript::calculate_helper_variables<std::vector<tobject>, tobject>(x0, aux);
		}
    };
}
#endif // cl_tape_examples_impl_cubic_spline_utils_hpp