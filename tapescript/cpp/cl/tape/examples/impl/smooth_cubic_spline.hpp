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

#ifndef cl_tape_examples_impl_cubic_spline_hpp
#define cl_tape_examples_impl_cubic_spline_hpp

#include <cl/tape/tape.hpp>
#include "impl/utils.hpp"
#include "impl/smooth_cubic_spline_utils.hpp"

namespace cl
{
    class smooth_cubic_spline
    {
    protected:

        int data_number;
        // Smoothing parameter.
        tobject lambda = 0.9;
        // Input points.
        std::vector<tobject> x_input;
        // Input values.
        std::vector<tobject> y_input;
        // Calculation util.
        smooth_cubic_spline_utils *spline_math = new smooth_cubic_spline_utils;
        // Initialization flag.
        bool is_initialized = false;
        // Calculation flag.
        bool is_calculated = false;

    public:

        smooth_cubic_spline(){}
        virtual ~smooth_cubic_spline(){}

        // Sets smoothing parameter.
        virtual void set_lambda(tobject l)
        {
            lambda = l;
        }

        // Returns input grid.
        virtual std::vector<tobject> get_x()
        {
            return x_input;
        }

		// Loads input data.
		virtual void load(const std::vector<double>& x, const std::vector<tobject>& y)
		{
			std::vector<tobject> tx;
			for each (auto e in x)
				tx.emplace_back(e);
			load(tx, y);
		}

        // Loads input data.
        virtual void load(const std::vector<tobject>& x, const std::vector<tobject>& y)
        {
            if (x.size() != y.size())
                throw std::exception("Wrong input data");
            x_input = x;
            y_input = y;
            data_number = y.size();
            (*spline_math).load(x, y);
            is_initialized = true;
        }

        // Calculates spline parameters.
        virtual void calculate()
		{
			if (!is_initialized)
				throw std::exception("Spline is not initialised.");
			// Estimate the optimal rigidity.
			//estimate_rigidity(); // for dependent lambda use set lambda after tape_start(X);
			// Solve underlying linear system to find spline coefficients.
			(*spline_math).Solve(lambda);
			is_calculated = true;
		}

		// Estimates the best rigidity value for smoothing spline.
		virtual void estimate_rigidity()
		{
			int n = y_input.size();
			// Rigidity value lies between 0 and 1.
			if (n < 1) lambda = 1e-6;
			std::vector<tobject> sT(n);
			// Estimating the variance of the set of delta y instead of y to consider the rotation of the function.
			for (int i = 0; i < n - 1; i++)
				sT[i] = y_input[i + 1] - y_input[i];
			tobject variance = 0;
			tobject t = sT[0];
			// Estimation of the variance.
			for (int i = 1; i < n; i++)
			{
				t += sT[i];
				tobject diff = ((i + 1) * sT[i]) - t;
				variance += (diff * diff) / ((i + 1.0) * i);
			}
			tobject varSt = variance / (n - 1);
			// Calculating the value of the rigidity.
			lambda = 1.0 / (varSt * std::pow(n, 1.0 / 2) + 1);
		}

        // Returns spline value at a specified point.
        virtual tobject value_at(tobject point)
        {
            if (!is_calculated)
                throw new std::exception("Spline is not calculated.");
            return (*spline_math).value_at(point);
        }

        virtual std::vector<tobject> fitted_values()
        {
            return (*spline_math).get_fitted_values();
        }

		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		// Returns lambda.
		virtual tobject get_lambda()
		{
			return lambda;
		}

		virtual std::vector<tobject> spline_vec(std::vector<tobject> x0, std::vector<tobject> y0)
		{
			load(x0, y0);
			calculate();
			return fitted_values();
		}
    };
}
#endif // cl_tape_examples_impl_cubic_spline_hpp