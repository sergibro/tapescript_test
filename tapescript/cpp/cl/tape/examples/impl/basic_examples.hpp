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

#ifndef cl_basic_examples_hpp
#define cl_basic_examples_hpp

#define CL_BASE_SERIALIZER_OPEN
#include <cl/tape/tape.hpp>

namespace cl
{
    inline void add_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Testing tape contents for addition of two tdoubles." << std::endl;

        // Input values initialization.
        std::vector<tdouble> X(2, 3.0);

        // Declare the X vector as independent and start tape recording.
        cl::tape_start(X);

        // Output calculations.
        std::vector<tdouble> Y(1);

        Y[0] = 2 * X[0] + X[1];

        out_str << "\nFunction Y = 2 * X[0] + X[1] is being tested at X[0] = " << X[0] << ", X[1] = " << X[1] << std::endl;

        // Declare a tape function and stop the tape recording.
        cl::tfunc<double> f(X, Y);

        // Calculate derivative of Y with respect to X[0] in Forward mode.
        std::vector<double> sy, sx(2, 0.0);
        sx[0] = 1;

        out_str << "\nTape operations sequence for differentiation of Y with respect to X[0]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[0] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << std::endl;

        // Calculate derivative of Y with respect to X[1] in Forward mode.
        sx[0] = 0;
        sx[1] = 1;

        out_str << "\nTape operations sequence for differentiation of Y with respect to X[1]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[1] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX(X[1]) = " << sy[0] << std::endl;

        // Calculate derivatives in Reverse mode.
        std::vector<double> sw(1, 1.0);

        out_str << "\nTape operations sequence for differentiation of Y with respect to X: ";

        sy = f.reverse(1, sw, out_stream);

        out_str << "\nDerivatives of Y with respect to X in Reverse mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << ", dY / dX (X[1]) = " << sy[1] << std::endl;
        out_str << std::endl << std::string(100, '-') << std::endl;
    }

    inline void mult_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Testing tape contents for multiplying of two tdoubles." << std::endl;

        // Input values initialization.
        std::vector<tdouble> X(2);
        X[0] = 2.0;
        X[1] = 3.0;

        // Declare the X vector as independent and start a tape recording.
        cl::tape_start(X);

        // Output calculations.
        std::vector<tdouble> Y(1);

        Y[0] = X[0] * X[1];

        out_str << "\nFunction Y = X[0] * X[1] is being tested at X[0] = " << X[0] << ", X[1] = " << X[1] << std::endl;

        // Declare a tape function and stop the tape recording.
        cl::tfunc<double> f(X, Y);

        // Calculate derivative of Y with respect to X[0] in Forward mode.
        std::vector<double> sy, sx(2, 0.0);
        sx[0] = 1;

        out_str << "\nTape operations sequence for differentiation of Y with respect to X[0]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[0] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX(X[0]) = " << sy[0] << std::endl;

        // Calculate derivative of Y on X[1] in Forward mode.
        sx[0] = 0;
        sx[1] = 1;

        out_str << "\nTape operations sequence for differentiation of Y with respect to X[1]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[1] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX(X[1]) = " << sy[0] << std::endl;

        // Calculate derivatives in Reverse mode.
        std::vector<double> sw(1, 1.0);

        out_str << "\nTape operations sequence for differentiation of Y with respect to X: ";

        sy = f.reverse(1, sw, out_stream);

        out_str << "\nDerivatives of Y on X in Reverse mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX(X[0]) = " << sy[0] << ", dY / dX(X[1]) = " << sy[1] << std::endl;
        out_str << std::endl << std::string(100, '-') << std::endl;
    }

    inline void pow_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Testing tape contents for pow of tdouble." << std::endl;

        // Input values initialization.
        std::vector<tdouble> X(1, 3.0);

        // Declare the X vector as independent and start a tape recording.
        cl::tape_start(X);

        // Output calculations.
        std::vector<tdouble> Y(1);

        Y[0] = std::pow(X[0], 4);

        out_str << "\nFunction Y =  X[0]^4 is being tested at X[0] = " << X[0] << std::endl;

        // Declare a tape function and stop the tape recording.
        cl::tfunc<double> f(X, Y);

        // Calculate derivative in Forward mode.
        std::vector<double> sy, sx(1, 1.0);

        out_str << "\nTape operations sequence for differentiation of Y with respect to X[0]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[0] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX(X[0]) = " << sy[0] << std::endl;

        // Calculate derivative in Reverse mode.
        std::vector<double> sw(1, 1.0);

        out_str << "\nTape operations sequence for differentiation of Y with respect to X: ";

        sy = f.reverse(1, sw, out_stream);

        out_str << "\nDerivative of Y with respect to X in Reverse mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << std::endl;
        out_str << std::endl << std::string(100, '-') << std::endl;
    }

    inline void exponent_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Testing tape contents for exponent of tdouble." << std::endl;

        // Input values initialization.
        std::vector<tdouble> X(1, 3.0);
        X[0] += std::exp(X[0]);
        tdouble x_exp = 3.0 + std::exp(3.0);
        CL_ASSERT(X[0] == x_exp, "Calculated and expected values are different.");

        // Declare the X vector as independent and start a tape recording.
        cl::tape_start(X);

        // Output calculations.
        std::vector<tdouble> Y(1);

        Y[0] = std::exp(-2.0 * X[0]);
        out_str << "\nFunction Y = exp(-2.0 * X[0]) is being tested at X[0] = " << X[0] << std::endl;

        // Declare a tape function and stop the tape recording.
        cl::tfunc<double> f(X, Y);

        // Calculate derivative in Forward mode.
        std::vector<double> sy, sx(1, 1.0);

        out_str << "\nTape operations sequence for differentiation of Y with respect to X[0]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[0] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << std::endl;

        // Calculate derivative in Reverse mode.
        std::vector<double> sw(1, 1.0);

        out_str << "\nTape operations sequence for differentiation of Y with respect to X: ";

        sy = f.reverse(1, sw, out_stream);

        out_str << "\nDerivative of Y with respect to X in Reverse mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << std::endl;
        out_str << std::endl << std::string(100, '-') << std::endl;
    }

    inline void log_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Testing tape contents for log of tdouble." << std::endl;

        // Input values initialization.
        std::vector<tdouble> X(2, 3.0);
        X[0] += std::log(X[0]);
        X[1] += std::log10(X[1]);
        tdouble x0_exp = 3.0 + std::log(3.0);
        tdouble x1_exp = 3.0 + std::log10(3.0);
        CL_ASSERT(X[0] == x0_exp && X[1] == x1_exp, "Calculated and expected values are different.");

        // Declare the X vector as independent and start a tape recording.
        cl::tape_start(X);

        // Output calculations.
        std::vector<tdouble> Y(2);

        Y[0] = std::log(2.0 * X[0]);
        Y[1] = std::log10(2.0 * X[1]);
        out_str << "\nFunction Y[0] = log(2.0 * X[0]) is being tested at X[0] = " << X[0] << std::endl;
        out_str << "\nFunction Y[1] = log10(2.0 * X[1]) is being tested at X[1] = " << X[1] << std::endl;

        // Declare a tape function and stop the tape recording.
        cl::tfunc<double> f(X, Y);

        // Calculate derivative in Forward mode.
        std::vector<double> sy, sx(2, 0.0);
        sx[0] = 1;
        out_str << "\nTape operations sequence for differentiation of Y with respect to X[0]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[0] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << std::endl;

        // Calculate derivative in Reverse mode.
        std::vector<double> sw(2, 0.0);
        sw[0] = 1;
        out_str << "\nTape operations sequence for differentiation of Y with respect to X: ";

        sy = f.reverse(1, sw, out_stream);

        out_str << "\nDerivative of Y with respect to X in Reverse mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << std::endl;
        out_str << std::endl << std::string(100, '-') << std::endl;
    }

    inline void abs_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Testing tape contents for abs of tdouble." << std::endl;

        // Input values initialization.
        std::vector<tdouble> X(1, -3.0);
        X[0] -= std::abs(X[0]);
        tdouble x_exp = -3.0 - std::abs(-3.0);
        CL_ASSERT(X[0] == x_exp, "Calculated and expected values are different.");

        // Declare the X vector as independent and start a tape recording.
        cl::tape_start(X);

        // Output calculations.
        std::vector<tdouble> Y(1);

        Y[0] = std::abs(2.0 * X[0]);
        out_str << "\nFunction Y = abs(2.0 * X[0]) is being tested at X[0] = " << X[0] << std::endl;

        // Declare a tape function and stop the tape recording.
        cl::tfunc<double> f(X, Y);

        // Calculate derivative in Forward mode.
        std::vector<double> sy, sx(1, 1.0);

        out_str << "\nTape operations sequence for differentiation of Y with respect to X[0]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[0] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << std::endl;

        // Calculate derivative in Reverse mode.
        std::vector<double> sw(1, 1.0);

        out_str << "\nTape operations sequence for differentiation of Y with respect to X: ";

        sy = f.reverse(1, sw, out_stream);

        out_str << "\nDerivative of Y with respect to X in Reverse mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << std::endl;
        out_str << std::endl << std::string(100, '-') << std::endl;
    }


    inline void trigonometric_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Testing tape contents for addition of cos and sin of tdouble." << std::endl;

        // Input values initialization.
        std::vector<tdouble> X(2);
        X[0] = 2.0;
        X[1] = 3.0;
        tdouble y_exp = std::cos(X[0]) + std::sin(X[1]) + std::tan(X[0] + X[1]);
        // Declare the X vector as independent and start a tape recording.
        cl::tape_start(X);

        // Output calculations.
        std::vector<tdouble> Y(1);

        Y[0] = std::cos(X[0]) + std::sin(X[1]) + std::tan(X[0]+ X[1]);

        out_str << "\nFunction Y = cos(X[0]) + sin(X[1] + tan(X[0]+ X[1])) is being tested at X[0] = " << X[0] << ", X[1] = " << X[1] << std::endl;

        // Declare a tape function and stop the tape recording.
        cl::tfunc<double> f(X, Y);

        // Calculate derivative of Y with respect to X[0] in Forward mode.
        std::vector<double> sy, sx(2, 0.0);
        sx[0] = 1;

        out_str << "\nTape operations sequence for differentiation of Y with respect to X[0]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[0] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << std::endl;

        // Calculate derivative of Y with respect to X[1] in Forward mode.
        sx[0] = 0;
        sx[1] = 1;

        out_str << "\nTape operations sequence for differentiation of Y with respect to X[1]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[1] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[1]) = " << sy[0] << std::endl;

        // Calculate derivatives in Reverse mode.
        std::vector<double> sw(1, 1.0);

        out_str << "\nTape operations sequence for differentiation of Y with respect to X: ";

        sy = f.reverse(1, sw, out_stream);

        out_str << "\nDerivative of Y with respect to X in Reverse mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << ", dY / dX (X[1]) = " << sy[1] << std::endl;
    }

    inline void hyperbolic_trigonometric_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Testing tape contents for addition of cosh, sinh and tanh of tdouble." << std::endl;

        // Input values initialization.
        std::vector<tdouble> X(2);
        X[0] = 2.0;
        X[1] = 3.0;
        tdouble y_exp = std::cosh(X[0]) + std::sinh(X[1]) + std::tanh(X[0] + X[1]);
        // Declare the X vector as independent and start a tape recording.
        cl::tape_start(X);

        // Output calculations.
        std::vector<tdouble> Y(1);

        Y[0] = std::cosh(X[0]) + std::sinh(X[1]) + std::tanh(X[0] + X[1]);

        out_str << "\nFunction Y = cosh(X[0]) + sinh(X[1] + tanh(X[0]+ X[1])) is being tested at X[0] = " << X[0] << ", X[1] = " << X[1] << std::endl;

        // Declare a tape function and stop the tape recording.
        cl::tfunc<double> f(X, Y);

        // Calculate derivative of Y with respect to X[0] in Forward mode.
        std::vector<double> sy, sx(2, 0.0);
        sx[0] = 1;

        out_str << "\nTape operations sequence for differentiation of Y with respect to X[0]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[0] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << std::endl;

        // Calculate derivative of Y with respect to X[1] in Forward mode.
        sx[0] = 0;
        sx[1] = 1;

        out_str << "\nTape operations sequence for differentiation of Y with respect to X[1]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[1] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[1]) = " << sy[0] << std::endl;

        // Calculate derivatives in Reverse mode.
        std::vector<double> sw(1, 1.0);

        out_str << "\nTape operations sequence for differentiation of Y with respect to X: ";

        sy = f.reverse(1, sw, out_stream);

        out_str << "\nDerivative of Y with respect to X in Reverse mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << ", dY / dX (X[1]) = " << sy[1] << std::endl;
    }

    inline void inverse_trigonometric_example(std::ostream& out_stream = std::cout)
    {
        out_str << "Testing tape contents for addition of acos, asin and atan of tdouble." << std::endl;

        // Input values initialization.
        std::vector<tdouble> X(2);
        X[0] = 0.3;
        X[1] = 0.2;
        tdouble y_exp = std::acos(X[0]) + std::asin(X[1]) + std::atan(X[0] + X[1]) + std::atan2(X[1],X[0]);
        // Declare the X vector as independent and start a tape recording.
        cl::tape_start(X);

        // Output calculations.
        std::vector<tdouble> Y(1);

        Y[0] = std::acos(X[0]) + std::asin(X[1]) + std::atan(X[0] + X[1]) + std::atan2(X[1], X[0]);

        out_str << "\nFunction Y = acos(X[0]) + asin(X[1] + atan(X[0]+ X[1]) + atan2(X[1],X[0])) is being tested at X[0] = " << X[0] << ", X[1] = " << X[1] << std::endl;

        // Declare a tape function and stop the tape recording.
        cl::tfunc<double> f(X, Y);

        // Calculate derivative of Y with respect to X[0] in Forward mode.
        std::vector<double> sy, sx(2, 0.0);
        sx[0] = 1;

        out_str << "\nTape operations sequence for differentiation of Y with respect to X[0]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[0] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << std::endl;

        // Calculate derivative of Y with respect to X[1] in Forward mode.
        sx[0] = 0;
        sx[1] = 1;

        out_str << "\nTape operations sequence for differentiation of Y with respect to X[1]: ";

        sy = f.forward(1, sx, out_stream);

        out_str << "Derivative of Y with respect to X[1] in Forward mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[1]) = " << sy[0] << std::endl;

        // Calculate derivatives in Reverse mode.
        std::vector<double> sw(1, 1.0);

        out_str << "\nTape operations sequence for differentiation of Y with respect to X: ";

        sy = f.reverse(1, sw, out_stream);

        out_str << "\nDerivative of Y with respect to X in Reverse mode has been calculated successfully:" << std::endl;
        out_str << "\tdY / dX (X[0]) = " << sy[0] << ", dY / dX (X[1]) = " << sy[1] << std::endl;
    }

    inline void basic_examples()
    {
        std::ofstream output("output/basic_examples_output.txt");
        cl::tape_serializer<double> serializer(output);

        add_example(serializer);
        mult_example(serializer);
        pow_example(serializer);
        exponent_example(serializer);
        log_example(serializer);
        abs_example(serializer);
        trigonometric_example(serializer);
        hyperbolic_trigonometric_example(serializer);
        inverse_trigonometric_example(serializer);

    }
}

#endif // cl_basic_examples_hpp
