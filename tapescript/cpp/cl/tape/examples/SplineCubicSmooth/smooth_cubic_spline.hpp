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
    private:

        int data_number;
        // Smoothing parameter.
        double lambda = 0.9;
        // Input points.
        std::vector<double> x_input;
        // Input values.
        std::vector<tdouble> y_input;
        // Calculation util.
        smooth_cubic_spline_utils spline_math;
        // Initialization flag.
        bool is_initialized = false;
        // Calculation flag.
        bool is_calculated = false;

    public:

        smooth_cubic_spline(){}
        ~smooth_cubic_spline(){}

        // Sets smoothing parameter.
        void set_lambda(double l)
        {
            lambda = l;
        }

        // Returns input grid.
        std::vector<double> get_x()
        {
            return x_input;
        }

        // Loads input data.
        void load(const std::vector<double>& x, const std::vector<tdouble>& y)
        {
            if (x.size() != y.size())
                throw std::exception("Wrong input data");
            x_input = x;
            y_input = y;
            data_number = y.size();
            spline_math.load(x, y);
            is_initialized = true;
        }

        // Calculates spline parameters.
        void calculate()
        {
            if (!is_initialized)
                throw std::exception("Spline is not initialised.");

            // Solve underlying linear system to find spline coefficients.
            spline_math.Solve(lambda);

            is_calculated = true;
        }

        // Returns spline value at a specified point.
        tdouble value_at(double point)
        {
            if (!is_calculated)
                throw new std::exception("Spline is not calculated.");
            return spline_math.value_at(point);
        }

        std::vector<tdouble> fitted_values()
        {
            return spline_math.get_fitted_values();
        }

    };
}
#endif // cl_tape_examples_impl_cubic_spline_hpp