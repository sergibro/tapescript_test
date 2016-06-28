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

#ifndef cl_tape_examples_impl_smooth_cubic_spline_example_hpp
#define cl_tape_examples_impl_smooth_cubic_spline_example_hpp

#include <cl/tape/tape.hpp>
#include "impl/utils.hpp"
#include "impl/smooth_cubic_spline.hpp"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

namespace cl
{
    tdouble func(tdouble x);

    tdouble deviation(smooth_cubic_spline* spline, std::vector<tdouble> y, std::function<tdouble(tdouble)> function);

    std::vector<std::vector<tdouble>> finite_difference_derivatives(const std::vector<double>& x, const std::vector<tdouble>& y, double step);

    void data(int n, std::vector<double>& x, std::vector<tdouble>& y, std::function<tdouble(tdouble)> function, boost::normal_distribution<double> nd);

    inline void smooth_cubic_spline_example(std::ostream& out_stream = std::cout)
    {
        const int n = 10;
        std::vector<double> x(n + 1);
        std::vector<tdouble> y(n + 1);
        boost::mt19937 rng;
        boost::normal_distribution<double> nd(0.0, 0);
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<>> var_nor(rng, nd);
        for (int i = 0; i < n + 1; i++)
        {
            x[i] = (double)i / n;
            tdouble noise = var_nor();
            y[i] = func(x[i]) + noise;
        }
        out_stream << "x " << x << "\n";
        out_stream << "y " << y << "\n";
        smooth_cubic_spline spline;
        spline.load(x, y);
        spline.calculate();
        std::vector<tdouble> y_f = spline.fitted_values();
        out_stream << "y_fitted " << y_f << "\n";
        tdouble dev = deviation(&spline, y_f, func);
        out_stream << std::endl << "deviation between estimation and reference: " << dev << std::endl << std::endl;
    }

    // Test for differentiation of smooth cubic spline.
    inline void smooth_cubic_spline_example_with_tdouble(std::ostream& out_stream = std::cout)
    {
        const int n = 9;
        double step = 1e-3;
        std::vector<double> x;
        std::vector<tdouble> y;
        boost::normal_distribution<double> nd(0.0, 5);
        data(n, x, y, func, nd);
        smooth_cubic_spline spline;
        tcontext<tdouble> spline_context;
        for (int i = 0; i < y.size(); i++)
        {
            std::stringstream ss;
            ss << i;
            std::string str("y_noised_");
            ss >> str;
            spline_context.independent(y[i], str);
        }
        std::vector<tdouble> y_context(n + 1);
        for (int i = 0; i < n + 1; i++)
            y_context[i] = spline_context[i];
        spline.load(x, y_context);
        spline.calculate();
        std::vector<tdouble> y_f = spline.fitted_values();
        for (int i = 0; i < y.size(); i++)
        {
            std::stringstream ss;
            ss << i;
            std::string str("y_fitted_");
            ss >> str;
            spline_context.dependent(y_f[i], str);
        }
        std::vector<double> forward;
        std::vector<double> reverse;
        std::vector<double> dx(y.size());
        long finite_difference_time = -std::clock();
        std::vector<std::vector<tdouble>> finite_difference = finite_difference_derivatives(x, y, step);
        finite_difference_time += std::clock();
        long forward_time = 0;
        long reverse_time = 0;
        out_stream << "Differentiating fitted y values by input y values... \n\n";
        for (int i = 0; i < y.size(); i++)
        {
            for (int j = 0; j < y.size(); j++)
                dx[j] = 0;
            dx[i] = 1;
            forward_time -= std::clock();
            forward = spline_context.forward(1, dx);
            forward_time += std::clock();
            reverse_time -= std::clock();
            reverse = spline_context.reverse(1, dx);
            reverse_time += std::clock();
            out_stream << "Forward derivatives = " << forward << "..." << std::endl;
            out_stream << "Reverse derivatives = " << reverse << "..." << std::endl;
            out_stream << "Finite-difference derivatives = " << finite_difference[i] << "..." << std::endl;
        }
        //out_stream << "Forward calculation took (ms): " << forward_time / (double)(CLOCKS_PER_SEC)* 1000 << std::endl;
        //out_stream << "Reverse calculation took (ms): " << reverse_time / (double)(CLOCKS_PER_SEC)* 1000 << std::endl;
        //out_stream << "Finite-difference calculation took (ms): " << finite_difference_time / (double)(CLOCKS_PER_SEC)* 1000 << std::endl;
    }

    // Calculates deviation between fitted values and reference.
    inline tdouble deviation(smooth_cubic_spline* spline, std::vector<tdouble> y, std::function<tdouble(tdouble)> function)
    {
        tdouble sum = 0;
        std::vector<double> x = spline->get_x();
        int n = x.size();
        std::vector<tdouble> y_fitted = spline->fitted_values();
        for (int i = 0; i < n; i++)
        {
            tdouble y_fitted_value = y[i];
            tdouble y_reference = function(x[i]);
            sum += (y_reference - y_fitted_value)*(y_reference - y_fitted_value);
        }
        //delete y_fitted;
        return std::sqrt(sum / n);
    }

    // Test function.
    inline tdouble func(tdouble x)
    {
        return x*x;
    }

    // Auxiliary method for generating noisy data.
    inline void data(int n, std::vector<double>& x, std::vector<tdouble>& y, std::function<tdouble(tdouble)> function, boost::normal_distribution<double> nd)
    {
        x.resize(n + 1);
        y.resize(n + 1);
        boost::mt19937 rng;
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<>> var_nor(rng, nd);
        for (int i = 0; i < n + 1; i++)
        {
            x[i] = (double)i / n;
            tdouble noise = var_nor();
            y[i] = function(x[i]) + noise;
        }
    }

    // Finite-difference smooth spline derivatives.
    std::vector<std::vector<tdouble>> finite_difference_derivatives(const std::vector<double>& x, const std::vector<tdouble>& y, double step)
    {
        int n = y.size();
        smooth_cubic_spline spline;
        spline.load(x, y);
        spline.calculate();
        std::vector<tdouble> first = spline.fitted_values();
        std::vector<tdouble> y_shifted(n);
        std::vector<std::vector<tdouble>> derivative(n);
        for (int i = 0; i < n; i++)
            y_shifted[i] = y[i];
        for (int i = 0; i < n; i++)
        {
            y_shifted[i] = y_shifted[i] + step;
            smooth_cubic_spline spline_next;
            spline_next.load(x, y_shifted);
            spline_next.calculate();
            std::vector<tdouble> second = spline_next.fitted_values();
            derivative[i].resize(n);
            for (int j = 0; j < n; j++)
            {
                tdouble diff = second[j] - first[j];
                derivative[i][j] = diff / step;
            }
            y_shifted[i] = y_shifted[i] - step;
        }
        return derivative;
    }

    inline void cubic_spline_examples()
    {
        std::ofstream of("output/smooth_cubic_spline_example_output.txt");
        cl::tape_serializer<cl::tvalue> serializer(of);
        serializer.precision(3);

        smooth_cubic_spline_example(serializer);
        smooth_cubic_spline_example_with_tdouble(serializer);
    }
}
#endif // cl_tape_examples_impl_smooth_cubic_spline_example_hpp