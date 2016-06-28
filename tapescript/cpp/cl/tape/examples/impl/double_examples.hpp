/*
Copyright (C) 2015-present CompatibL

Performance test results and finance-specific examples are available at:

http://www.tapescript.org

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

#ifndef cl_tape_examples_impl_double_examples_hpp
#define cl_tape_examples_impl_double_examples_hpp

#define CL_BASE_SERIALIZER_OPEN
#include <cl/tape/tape.hpp>
#include "impl/array_examples.hpp"


namespace cl
{
    inline void plus_double_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Plus operation:\n" << std::endl;

        // Initialize input values
        std::vector<tdouble> X = { 1, 2, 0, 1 };
        out_str << "Input vector: " << X << "\n";

        // Declare the X vector as independent and start a tape recording.
        tape_start(X);

        // Output calculations.
        std::vector<tdouble> sum = { X[0] + X[2], X[1] + X[3] };
        std::vector<tdouble> Y = sum;
        out_str << "Output vector: " << Y << "\n\n";

        // Declare a tape function and stop the tape recording.
        tfunc<double> f(X, Y);

        // Forward sweep calculations.
        std::vector<double> dx = { 2, 5, 1, -1 };
        out_str << "Forward(1, dx) sweep for dx = " << dx << "..." << std::endl;
        std::vector<double> forw = f.forward(1, dx, out_stream);
        out_str << "Forward sweep result: " << forw << "\n\n";

        // Reverse sweep calculations.
        std::vector<double> w =  { 3, 7 } ;
        out_str << "Reverse(1, w) sweep for w = " << w << "..." << std::endl;
        std::vector<double> rev = f.reverse(1, w);
        out_str << "Reverse sweep result: " << rev << "\n\n\n";
    }

    inline void minus_double_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Minus operation:\n\n";

        // Input values initialization.
        std::vector<tdouble> X = { 1, 2, 0, 1 };
        out_str << "Input vector: " << X << "\n";

        // Declare the X vector as independent and start a tape recording.
        tape_start(X);

        // Output calculations.
        std::vector<tdouble> diff = { X[0] - X[2], X[1] - X[3] };
        std::vector<tdouble> Y = diff;
        out_str << "Output vector: " << Y << "\n\n";

        // Declare a tape function and stop the tape recording.
        tfunc<double> f(X, Y);

        // Forward sweep calculations.
        std::vector<double> dx = { 2, 5, 1, -1 };
        out_str << "Forward(1, dx) sweep for dx = " << dx << "..." << std::endl;
        std::vector<double> forw = f.forward(1, dx, out_stream);
        out_str << "Forward sweep result: " << forw << "\n\n";

        // Reverse sweep calculations.
        std::vector<double> w = { 3, 7 };
        out_str << "Reverse(1, w) sweep for w = " << w << "..." << std::endl;
        std::vector<double> rev = f.reverse(1, w);
        out_str << "Reverse sweep result: " << rev << "\n\n\n";
    }

    inline void exp_double_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Exponent function:\n\n";

        // Initialize input values
        std::vector<tdouble> X = { 0, 1 };
        out_str << "Input vector: " << X << "\n";

        // Declare the X vector as independent and start a tape recording.
        tape_start(X);

        // Output calculations.
        std::vector<tdouble> exp = { std::exp(X[0]), std::exp(X[1]) };
        std::vector<tdouble> Y = exp;
        out_str << "Output vector: " << Y << "\n\n";

        out_str << "Initial Forward(0) sweep...\n\n";
        // Declare a tape function and stop the tape recording.
        tfunc<double> f(X, Y);

        // Forward sweep calculations.
        std::vector<double> dx = { { 2, -1 } };
        out_str << "Forward(1, dx) sweep for dx = " << dx << "..." << std::endl;
        std::vector<double> forw = f.forward(1, dx, out_stream);
        out_str << "Forward sweep result: " << forw << "\n\n";

        // Reverse sweep calculations.
        std::vector<double> w = { { 1, -1 } };
        out_str << "Reverse(1, w) sweep for w = " << w << "..." << std::endl;
        std::vector<double> rev = f.reverse(1, w);
        out_str << "Reverse sweep result: " << rev << "\n\n\n";
    }

    inline void cos_double_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Cosine function:\n\n";

        // Initialize input values
        std::vector<tdouble> X = { 1, 3.14159265359 / 2 };
        out_str << "Input vector: " << X << "\n";

        // Declare the X vector as independent and start a tape recording.
        tape_start(X);

        // Output calculations.
        std::vector<tdouble> cos = { std::cos(X[0]), std::cos(X[1]) };
        std::vector<tdouble> Y = cos;
        out_str << "Output vector: " << Y << "\n\n";

        out_str << "Initial Forward(0) sweep...\n\n";
        // Declare a tape function and stop the tape recording.
        tfunc<double> f(X, Y);

        // Forward sweep calculations.
        std::vector<double> dx = { { 1, 1 } };
        out_str << "Forward(1, dx) sweep for dx = " << dx << "..." << std::endl;
        std::vector<double> forw = f.forward(1, dx, out_stream);
        out_str << "Forward sweep result: " << forw << "\n\n";

        // Reverse sweep calculations.
        std::vector<double> w = { { -2, 1 } };
        out_str << "Reverse(1, w) sweep for w = " << w << "..." << std::endl;
        std::vector<double> rev = f.reverse(1, w);
        out_str << "Reverse sweep result: " << rev << "\n\n\n";
    }

    inline void discount_double_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Discount calculations:\n\n";

        // Input values initialization.
        std::vector<tdouble>  spot = { 100, 110 };
        std::vector<tdouble>  rate = { 0.05, 0.04 };
        std::vector<tdouble>  time = { 1, 2 };
        std::vector<tdouble> X = { spot[0], spot[1], rate[0], rate[1], time[0], time[1] };
        out_str << "Spot: " << spot << "\n";
        out_str << "Rate: " << rate << "\n";
        out_str << "Time: " << time << "\n";
        out_str << "Input vector: " << X << "\n";

        // Declare the X vector as independent and start a tape recording.
        tape_start(X);

        // Output calculations.
        std::vector<tdouble> discount_factor = { std::exp(-X[2] * X[4]), std::exp(-X[3] * X[5]) };
        std::vector<tdouble> discounted_price = { X[0] * discount_factor[0], X[1] * discount_factor[1] };
        std::vector<tdouble> Y = { discounted_price[0], discounted_price[1], discount_factor[0], discount_factor[1] };
        out_str << "Discounted price: " << discounted_price << "\n";
        out_str << "Discount factor: " << discount_factor << "\n";
        out_str << "Output vector: " << Y << "\n\n";

        out_str << "Initial Forward(0) sweep...\n\n";
        // Declare a tape function and stop the tape recording.
        tfunc<double> f(X, Y);

        // Forward sweep calculations.
        out_str << "Calculation of sensitivity to 1% change in interest rate.\n";
        std::vector<double> dx = { 0, 0, 0.01, 0.01, 0, 0 };
        out_str << "Forward(1, dx) sweep for dx = " << dx << "..." << std::endl;
        std::vector<double> forw = f.forward(1, dx, out_stream);
        out_str << "Forward sweep result: " << forw << "\n";
        out_str << "Sensitivity of discounted price: " << forw[0] << " , " << forw[1] << "\n";
        out_str << "Sensitivity of discount factor:  " << forw[2] << " , " << forw[3] << "\n\n";

        // Reverse sweep calculations.
        out_str << "Calculation of sensitivity of discount factor.\n";
        std::vector<double> w = { 0, 0, 1, 1 };
        out_str << "Reverse(1, w) sweep for w = " << w << "..." << std::endl;
        std::vector<double> rev = f.reverse(1, w);
        out_str << "Reverse sweep result: " << rev << "\n";
        out_str << "Sensitivity to the changes in spot price:    " << rev[0] << " , " << rev[1] << "\n";
        out_str << "Sensitivity to the changes in interest rate: " << rev[2] << " , " << rev[3] << "\n";
        out_str << "Sensitivity to the changes in maturity time: " << rev[4] << " , " << rev[5] << "\n\n\n";
    }

    inline void double_example1(std::ostream& out_stream = std::cout)
    {
        out_str << "We are mixing arrays and scalars in calculations of\n";
        out_str << "Y = { sin(x0) / cos(x1) + x1 / x0, y0 + x0 * x1^3 - 0.5 } function:\n\n";

        // Initialize input values
        std::vector<tdouble> X = { 1, 2, -1.2 };
        out_str << "Input vector: " << X << "\n";

        // Declare the X vector as independent and start a tape recording.
        tape_start(X);

        // Output calculations.
        std::vector<tdouble> y0 = { std::sin(X[0]) / std::cos(X[2]) + X[2] / X[0], std::sin(X[1]) / std::cos(X[2]) + X[2] / X[1] };
        std::vector<tdouble> y1 = { y0[0] + X[0] * std::pow(X[2], 3) - 0.5, y0[1] + X[1] * std::pow(X[2], 3) - 0.5 };
        std::vector<tdouble> Y = { y0[0], y0[1], y1[0], y1[1] };
        out_str << "Output vector: " << Y << "\n\n";

        out_str << "Initial Forward(0) sweep...\n\n";
        // Declare a tape function and stop the tape recording.
        tfunc<double> f(X, Y);

        // Forward sweep calculations.
        std::vector<double> dx = { 1, -1 , 2 };
        out_str << "Forward(1, dx) sweep for dx = " << dx << "..." << std::endl;
        std::vector<double> forw = f.forward(1, dx, out_stream);
        out_str << "Forward sweep result: " << forw << "\n\n";

        // Reverse sweep calculations.
        std::vector<double> w = {  -2, 1, 2, 3 };
        out_str << "Reverse(1, w) sweep for w = " << w << "..." << std::endl;
        std::vector<double> rev = f.reverse(1, w);
        out_str << "Reverse sweep result: " << rev << "\n\n\n";
    }

    inline void sum_double_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Sum of an array:\n\n";

        // Initialize input values
        std::vector<tdouble> X = { 1, 2 };
        out_str << "Input vector: " << X << "\n";

        // Declare the X vector as independent and start a tape recording.
        tape_start(X);

        // Output calculations.
        tdouble y0 = X[0] + X[1];
        std::vector<tdouble> Y = { y0 };
        out_str << "Output vector: " << Y << "\n\n";

        out_str << "Initial Forward(0) sweep...\n\n";
        // Declare a tape function and stop the tape recording.
        tfunc<double> f(X, Y);

        // Forward sweep calculations.
        std::vector<double> dx = { 1.5, -1 };
        out_str << "Forward(1, dx) sweep for dx = " << dx << "..." << std::endl;
        std::vector<double> forw = f.forward(1, dx, out_stream);
        out_str << "Forward sweep result: " << forw << "\n\n";

        // Reverse sweep calculations.
        std::vector<double> w = { 4 };
        out_str << "Reverse(1, w) sweep for w = " << w << "..." << std::endl;
        std::vector<double> rev = f.reverse(1, w);
        out_str << "Reverse sweep result: " << rev << "\n\n\n";
    }

    inline void rev_double_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Reversing of an array:\n\n";

        // Input values initialization.
        std::vector<tdouble> X = { 1, 2 };
        out_str << "Input vector: " << X << "\n";

        // Declare the X vector as independent and start a tape recording.
        tape_start(X);

        // Output calculations.
        std::vector<tdouble> y0 = { X[1], X[0] };
        std::vector<tdouble > Y =  y0;
        out_str << "Output vector: " << Y << "\n\n";

        out_str << "Initial Forward(0) sweep...\n\n";
        // Declare a tape function and stop the tape recording.
        tfunc<double> f(X, Y);

        // Forward sweep calculations.
        std::vector<double> dx = { 1.5, -1 };
        out_str << "Forward(1, dx) sweep for dx = " << dx << "..." << std::endl;
        std::vector<double> forw = f.forward(1, dx, out_stream);
        out_str << "Forward sweep result: " << forw << "\n\n";

        // Reverse sweep calculations.
        std::vector<double> w = { 4, 4 };
        out_str << "Reverse(1, w) sweep for w = " << w << "..." << std::endl;
        std::vector<double> rev = f.reverse(1, w);
        out_str << "Reverse sweep result: " << rev << "\n\n\n";
    }

    inline void conc_double_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Array concatenation:\n" << std::endl;

        // Input values initialization.
        std::vector<tdouble> X = { 1, 2, 7, 0, 5 };
        out_str << "Input vector: " << X << "\n";

        // Declare the X vector as independent and start a tape recording.
        tape_start(X);

        // Calculate output
        std::vector<tdouble> Y = X;
        out_str << "Output vector: " << Y << "\n\n";

        out_str << "Initial Forward(0) sweep...\n\n";
        // Declare a tape function and stop the tape recording.
        tfunc<double> f(X, Y);

        // Forward sweep calculations.
        std::vector<double> dx = { 2, 5, 12, 1, -1 };
        out_str << "Forward(1, dx) sweep for dx = " << dx << "..." << std::endl;
        std::vector<double> forw = f.forward(1, dx, out_stream);
        out_str << "Forward sweep result: " << forw << "\n\n";

        // Reverse sweep calculations.
        std::vector<double> w = { 1, 2, 3, 4, 5 };
        out_str << "Reverse(1, w) sweep for w = " << w << "..." << std::endl;
        std::vector<double> rev = f.reverse(1, w);
        out_str << "Reverse sweep result: " << rev << "\n\n\n";
    }

    inline void make_double_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Array construction from scalar:\n" << std::endl;

        // Initialize input values
        std::vector<tdouble> X = { 7 };
        out_str << "Input vector: " << X << "\n";

        // Declare the X vector as independent and start a tape recording.
        tape_start(X);

        // Output calculations.
        std::vector<tdouble> y0 = { X[0], X[0] };
        std::vector<tdouble> Y = y0;
        out_str << "Output vector: " << Y << "\n\n";

        out_str << "Initial Forward(0) sweep...\n\n";
        // Declare a tape function and stop the tape recording.
        tfunc<double> f(X, Y);

        // Forward sweep calculations.
        std::vector<double> dx = { 1 };
        out_str << "Forward(1, dx) sweep for dx = " << dx << "..." << std::endl;
        std::vector<double> forw = f.forward(1, dx, out_stream);
        out_str << "Forward sweep result: " << forw << "\n\n";

        // Reverse sweep calculations.
        std::vector<double> w = { { 1, 2 } };
        out_str << "Reverse(1, w) sweep for w = " << w << "..." << std::endl;
        std::vector<double> rev = f.reverse(1, w);
        out_str << "Reverse sweep result: " << rev << "\n\n\n";
    }

    inline void max_double_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Maximum of two arrays:\n" << std::endl;

        // Initialize input values
        std::vector<tdouble> X = { -1, 0, 1, 1, 0, -1 };
        out_str << "Input vector: " << X << "\n";

        // Declare the X vector as independent and start a tape recording.
        tape_start(X);

        // Output calculations.
        std::vector<tdouble> y0 = { std::max(X[0], X[3]), std::max(X[1], X[4]), std::max(X[2], X[5]) };
        std::vector<tdouble> Y = y0;
        out_str << "Output vector: " << Y << "\n\n";

        out_str << "Initial Forward(0) sweep...\n\n";
        // Declare a tape function and stop the tape recording.
        tfunc<double> f(X, Y);

        // Forward sweep calculations.
        std::vector<double> dx = { 1, 1, 1, 2, 2, 2 };
        out_str << "Forward(1, dx) sweep for dx = " << dx << "..." << std::endl;
        std::vector<double> forw = f.forward(1, dx, out_stream);
        out_str << "Forward sweep result: " << forw << "\n\n";

        // Reverse sweep calculations.
        std::vector<double> w = { 1, 2, 3 };
        out_str << "Reverse(1, w) sweep for w = " << w << "..." << std::endl;
        std::vector<double> rev = f.reverse(1, w);
        out_str << "Reverse sweep result: " << rev << "\n\n\n";
    }

    inline void norm_double_example(std::ostream& out_stream = std::cout)
    {
        out_str << "L2 norm of an array:\n" << std::endl;

        // Initialize input values
        std::vector<tdouble> X = { 3 };
        out_str << "Input vector: " << X << "\n";

        // Declare the X vector as independent and start a tape recording.
        tape_start(X);

        // Output calculations.
        std::vector<tdouble> array = { X[0], X[0] };
        tdouble norm = std::sqrt(std::pow(array[0], 2) + std::pow(array[1], 2));
        std::vector<tdouble> Y = { norm };
        out_str << "Output vector: " << Y << "\n\n";

        out_str << "Initial Forward(0) sweep...\n\n";
        // Declare a tape function and stop the tape recording.
        tfunc<double> f(X, Y);

        // Forward sweep calculations.
        std::vector<double> dx = { 1 };
        out_str << "Forward(1, dx) sweep for dx = " << dx << "..." << std::endl;
        std::vector<double> forw = f.forward(1, dx, out_stream);
        out_str << "Forward sweep result: " << forw << "\n\n";

        // Reverse sweep calculations.
        std::vector<double> w = { 1 };
        out_str << "Reverse(1, w) sweep for w = " << w << "..." << std::endl;
        std::vector<double> rev = f.reverse(1, w);
        out_str << "Reverse sweep result: " << rev << "\n\n\n";
    }

    inline void linear_regression_double_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Linear regression:\n" << std::endl;

        // Input values initialization.
        size_t n = 3;
        std::vector<double> x_ref = { -1, 0, 1 };
        std::vector<double> y_ref = { 1, 0, 3 };
        std::vector<tdouble> X(2 * n);
        std::copy(x_ref.begin(), x_ref.end(), X.begin());
        std::copy(y_ref.begin(), y_ref.end(), X.begin() + n);

        out_str << "Input vector: " << X << "\n";

        // Declare the X vector as independent and start a tape recording.
        tape_start(X);

        // Output calculations.
        auto x = X.begin();
        auto y = X.begin() + n;

        tdouble x_mean = 0;
        tdouble y_mean = 0;
        for (size_t i = 0; i < n; i++)
        {
            x_mean += x[i];
            y_mean += y[i];
        }
        x_mean /= n;
        y_mean /= n;

        // Variance times n: n * Var[x]
        tdouble var_x_n = 0;
        // Covariance times n: n * Cov[x, y]
        tdouble cov_xy_n = 0;
        for (size_t i = 0; i < n; i++)
        {
            var_x_n += std::pow(x[i] - x_mean, 2);
            cov_xy_n += (x[i] - x_mean) * (y[i] - y_mean);
        }

        tdouble beta = cov_xy_n / var_x_n;
        tdouble alpha = y_mean - beta * x_mean;
        std::vector<tdouble> y_estimate(n);
        for (size_t i = 0; i < n; i++)
        {
            y_estimate[i] = alpha + beta * x[i];
        }
        std::vector<tdouble> Y(n + 2);
        Y[0] = alpha;
        Y[1] = beta;
        for (size_t i = 0; i < n; i++)
        {
            Y[i + 2] = y_estimate[i];
        }
        out_str << "Output vector: " << Y << "\n\n";

        out_str << "Initial Forward(0) sweep...\n\n";
        // Declare a tape function and stop the tape recording.
        tfunc<double> f(X, Y);

        // Forward sweep calculations.
        std::vector<double> dx = { 0, 0, 0, 0, 0, 1 };
        out_str << "Forward(1, dx) sweep for dx = " << dx << "..." << std::endl;
        std::vector<double> forw = f.forward(1, dx, out_stream);
        out_str << "Forward sweep result: " << forw << "\n\n";
    }

    inline void non_optimized_array_examples()
    {
        std::ofstream of("output/output_non_optimized_array.txt");
        tape_serializer<tvalue> serializer(of);
        serializer.precision(3);

        plus_double_example(serializer);
        minus_double_example(serializer);
        exp_double_example(serializer);
        cos_double_example(serializer);
        discount_double_example(serializer);
        double_example1(serializer);
        sum_double_example(serializer);
        rev_double_example(serializer);
        conc_double_example(serializer);
        make_double_example(serializer);
        max_double_example(serializer);
        norm_double_example(serializer);
        linear_regression_double_example(serializer);
    }
}

#endif // cl_tape_examples_impl_double_examples_hpp
