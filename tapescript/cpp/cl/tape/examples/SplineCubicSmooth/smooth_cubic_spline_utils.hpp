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
    private:

        // Number of input points.
        int data_number;
        // Spline sections.
        std::vector<double> x;
        std::vector<tdouble> y;
        // Auxiliary variables.
        std::vector<std::vector<tdouble>> auxiliaries;
        // Spline coefficients.
        std::vector<std::vector<tdouble>> coefficients;
        // Fitted values.
        std::vector<tdouble> y_fitted;

    public:

        smooth_cubic_spline_utils(){}
        ~smooth_cubic_spline_utils(){}

        std::vector<tdouble> get_fitted_values()
        {
            return y_fitted;
        }
        // Spline initialisation from X, Y input arrays.
        void load(const std::vector<double>& x0, const std::vector<tdouble>& y0)
        {
            // Number of points.
            data_number = x0.size();

            // Initialise spline input data.
            x = x0;
            y = y0;
        }

        // Calculate spline coefficients.
        void Solve(double lambda)
        {
            coefficients.resize(3);
            for (int i = 0; i < 3; i++)
                coefficients[i].resize(data_number - 1);
            // Calculate spline helper variables.
            cl::tapescript::calculate_helper_variables<std::vector<double>, tdouble>(x, auxiliaries);
            cl::tapescript::calculate_q<std::vector<tdouble>, tdouble>(y, auxiliaries);
            // Bound smoothing parameter value to avoid problems with NaN.
            double lambdaMin = 1e-7;
            double lambdaMax = 1 - 1e-7;
            if (lambda < lambdaMin)
                lambda = lambdaMin;
            else if (lambda > lambdaMax)
                lambda = lambdaMax;

            // Internal smoothing parameter: transform lambda -> mu = 2 * (1 - lambda) / (3 * lambda).
            double mu = 2 * (1 - lambda) / (3 * lambda);

            // Prepare system of linear equations: calculate diagonal, first and second upper diagonal vectors u, v, w, respectively.
            int n = data_number;
            std::vector<std::vector<tdouble>> uvw_;
            cl::tapescript::prepare_linear_system<tdouble>(mu, uvw_, auxiliaries);

            // Find b-coefficients in splines in auxiliaries[4] array.
            cl::tapescript::find_b_coefficients<tdouble>(uvw_, auxiliaries);
            std::vector<tdouble> d_coeffs = cl::tapescript::find_d_coefficients<std::vector<tdouble>, tdouble>(mu, y, auxiliaries);
            // Find the rest of spline parameters.
            cl::tapescript::find_all_coefficients<tdouble>(mu, d_coeffs, auxiliaries, coefficients);
            coefficients.push_back(d_coeffs);
            fitted_values();
        }

        void fitted_values()
        {
            int n = data_number;
            y_fitted.resize(n);
            for (int i = 0; i < n; i++)
                y_fitted[i] = value_at(x[i]);
        }

        /// <summary>Get spline derivative at point x with possible extrapolation.</summary>
        tdouble value_at(double point)
        {
            int n = data_number;
            if (point < x[0])
            {
                double diff = point - x[0];
                return coefficients[0][0] * std::pow(diff, 3) + coefficients[1][0] * std::pow(diff, 2) +
                    coefficients[2][0] * diff + coefficients[3][0];
            }
            else if (point > x[n - 1])
            {
                double diff = point - x[n - 1];
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

        tdouble value_in_section(double point, int section)
        {
            double diff = point - x[section];
            return coefficients[0][section] * std::pow(diff, 3) + coefficients[1][section] * std::pow(diff, 2) +
                coefficients[2][section] * diff + coefficients[3][section];
        }
    };
}
#endif // cl_tape_examples_impl_cubic_spline_utils_hpp