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

#include "impl/spline_cubic_smooth_end_value.hpp"
#include "impl/spline_cubic_smooth_slope.hpp"
#include "impl/spline_cubic_smooth_slope_value.hpp"

#define eps_test 1e-5
#define diff_with_respect_to_y  true

namespace cl
{
	std::vector<std::string> spline_test_data_file_names{ "11_put_geom_brownian", "13_put_geom_brownian", "15_put_geom_brownian", "17_straddle_geom_brownian", "19_straddle_geom_brownian",
		"1_call_geom_brownianGab", "1_RiskReversal", "1_Straddle", "1_Strangle", "21_straddle_geom_brownian",
		"23_straddle_geom_brownian", "25_risk reversal_geom_brownian", "27_risk reversal_geom_brownian", "29_risk reversal_geom_brownian", "2_RiskReversal",
		"2_Straddle", "2_Strangle", "31_risk reversal_geom_brownian", "3_put_geom_brownianGab", "3_RiskReversal",
		"3_Straddle", "3_Strangle", "4_RiskReversal", "4_Straddle", "4_Strangle",
		"5_RiskReversal", "7_call_geom_brownian", "9_put_geom_brownian", "AMCDataBear1D1000", "AMCDataBear1D10000",
		"AMCDataBear1D4000", "AMCDataBear1D7000", "AMCDataBull1D1000", "AMCDataBull1D10000", "AMCDataBull1D4000",
		"AMCDataBull1D7000", "AMCDataButterfly1D1000", "AMCDataButterfly1D10000", "AMCDataButterfly1D4000", "AMCDataButterfly1D7000",
		"DownAndInCallBarrierTest1000", "DownAndInCallBarrierTest10000", "DownAndOutCallBarrierTest1000", "DownAndOutCallBarrierTest10000", "PayoffCallOptionSell",
		"UpAndInCallBarrierTest1000", "UpAndInCallBarrierTest10000", "UpAndOutCallBarrierTest1000", "UpAndOutCallBarrierTest10000" };

	template<class T>
	std::string convertToStr(T *var)
	{
		std::ostringstream ss;
		ss << *var;
		return ss.str();
	}

	// Reads data from .csv file.
	template <class T>
	inline void read_file(std::string name, std::vector<std::vector<T>>& x, std::vector<T>& y, std::vector<T>& ref)
	{
		std::ifstream in("input/" + name + ".csv");
		std::string a;
		int n = 0;
		while (!in.eof())
		{
			in >> a;
			n++;
		}
		n--;
		int dim = 0;
		for (int i = 0; i < a.size(); i++)
		if (a[i] == 44)
			dim++;
		--dim;
		x.resize(n);
		for (int i = 0; i < n; i++)
			x[i].resize(dim);
		y.resize(n);
		ref.resize(n);
		in.close();
		std::ifstream in1("input/" + name + ".csv");
		for (int i = 0; i < n; i++)
		{
			in1 >> a;
			if (!i && (a[a[0] == '-' ? 1 : 0] > '9' || a[a[0] == '-' ? 1 : 0] < '0')) std::cout << std::endl << "Not a number! " + a << std::endl << std::endl; // check for correct read data
			std::string value = "";
			int d = 0;
			for (int j = 0; j < a.size(); j++)
			{
				if (a[j] == 44)
				{
					if (d < dim) x[i][d++] = (T)atof(value.c_str());
					else y[i] = (T)atof(value.c_str());
					value = "";
				}
				else value += a[j];
				if (j == (a.size() - 1)) ref[i] = (T)atof(value.c_str());
			}
		}
		in1.close();
	}
	
	template <class T>
	inline void read_and_prepare_data(std::string name, std::vector<T>& x, std::vector<T>& y, std::vector<T>& ref)
	{
		std::vector<std::vector<T>> xx;
		// Load data from file.
		read_file(name, xx, y, ref);

		for (auto e : xx) x.push_back(e.at(0));

		std::vector<std::vector<T>> for_sort(ref.size());
		for (auto j = 0; j < ref.size(); ++j)
		{
			for_sort.at(j).push_back(x.at(j));
			for_sort.at(j).push_back(y.at(j));
			for_sort.at(j).push_back(ref.at(j));
		}
		sort(for_sort.begin(), for_sort.end());
		x.clear();
		y.clear();
		ref.clear();
		auto x_tmp = for_sort.at(0).at(0);
		for (int j = 0; j < for_sort.size(); ++j)
		{
			if (for_sort.at(j).at(0) - x_tmp > eps_test || !j)
			{
				x.push_back(for_sort.at(j).at(0));
				y.push_back(for_sort.at(j).at(1));
				ref.push_back(for_sort.at(j).at(2));
				x_tmp = for_sort.at(j).at(0);
			}
		}
	}
	
	template <class T>
	void print_csv(std::vector<T> const& x, std::vector<T> const& y, std::vector<T> const& ref, std::vector<T> const& y_res, std::ostream& out = std::cout)
	{
		for (auto i = 1; i < y_res.size(); ++i)
			out << x.at(i) << ',' << y.at(i) << ',' << ref.at(i) << ',' << y_res.at(i) << std::endl;
	}
	
	template <class T>
	T delta(const std::vector<T>& y_res, const std::vector<T>& ref)
	{
		auto n = ref.size();
		T d = 0;
		for (auto i = 0; i < n; ++i)
		{
			d += (y_res.at(i) - ref.at(i)) * (y_res.at(i) - ref.at(i));
		}
		return std::sqrt(d / n);
	}

	// Estimates the best rigidity value for smoothing spline.
	template <class T>
	T estimate_rigidity(const std::vector<T>& y)
	{
		T lambda = 0.9;
		int n = y.size();
		// Rigidity value lies between 0 and 1.
		if (n < 1) lambda = 1e-6;
		std::vector<T> sT(n);
		// Estimating the variance of the set of delta y instead of y to consider the rotation of the function.
		for (int i = 0; i < n - 1; i++)
			sT[i] = y[i + 1] - y[i];
		T variance = 0;
		T t = sT[0];
		// Estimation of the variance.
		for (int i = 1; i < n; i++)
		{
			t += sT[i];
			T diff = ((i + 1) * sT[i]) - t;
			variance += (diff * diff) / ((i + 1.0) * i);
		}
		T varSt = variance / (n - 1);
		// Calculating the value of the rigidity.
		lambda = 1.0 / (varSt * std::pow(n, 1.0 / 2) + 1);
		return lambda;
	}

	// Test function.
	inline tobject func(tobject x)
	{
		return x*x;
	}

	// Calculates deviation between fitted values and reference.
	inline tobject deviation(smooth_cubic_spline* spline, std::vector<tobject> y, std::function<tobject(tobject)> function)
	{
		tobject sum = 0;
		std::vector<tobject> x = spline->get_x();
		int n = x.size();
		std::vector<tobject> y_fitted = spline->fitted_values();
		for (int i = 0; i < n; i++)
		{
			tobject y_fitted_value = y[i];
			tobject y_reference = function(x[i]);
			sum += (y_reference - y_fitted_value)*(y_reference - y_fitted_value);
		}
		//delete y_fitted;
		return std::sqrt(sum / n);
	}

	// Finite-difference smooth spline derivatives.
	std::vector<std::vector<tobject>> finite_difference_derivatives(const std::vector<double>& x, const std::vector<tobject>& y, double step)
	{
		int n = y.size();
		smooth_cubic_spline spline;
		spline.load(x, y);
		spline.calculate();
		std::vector<tobject> first = spline.fitted_values();
		std::vector<tobject> y_shifted(n);
		std::vector<std::vector<tobject>> derivative(n);
		for (int i = 0; i < n; i++)
			y_shifted[i] = y[i];
		for (int i = 0; i < n; i++)
		{
			y_shifted[i] = y_shifted[i] + step;
			smooth_cubic_spline spline_next;
			spline_next.load(x, y_shifted);
			spline_next.calculate();
			std::vector<tobject> second = spline_next.fitted_values();
			derivative[i].resize(n);
			for (int j = 0; j < n; j++)
			{
				tobject diff = second[j] - first[j];
				derivative[i][j] = diff / step;
			}
			y_shifted[i] = y_shifted[i] - step;
		}
		return derivative;
	}

	// Finite-difference smooth spline derivatives.
	template<typename T, typename TSpline>
	std::vector<std::vector<T>> finite_difference_derivatives(TSpline spline, std::vector<T>& x, std::vector<T> y, const int& point_number, const double& step = 1e-3)
	{
		std::vector<std::vector<T>> derivative(point_number);
		int n = y.size();
		int n_mult = n / (point_number + 1); // if end or slope value here we can't change 0 and n - 1 points
		for (int i = 0; i < point_number; ++i)
		{
			if (diff_with_respect_to_y) y[(i + 1) * n_mult] += step / 2;
			else x[(i + 1) * n_mult] += step / 2;
			std::vector<T> right = spline.spline_vec(x, y);
			if (diff_with_respect_to_y) y[(i + 1) * n_mult] -= step;
			else x[(i + 1) * n_mult] -= step;
			std::vector<T> left = spline.spline_vec(x, y);
			if (diff_with_respect_to_y) y[(i + 1) * n_mult] += step / 2;
			else x[(i + 1) * n_mult] += step / 2;
			derivative[i].resize(point_number);
			for (int j = 0; j < point_number; j++)
				derivative[i][j] = (right[(j + 1) * n_mult] - left[(j + 1) * n_mult]) / step;
		}
		return derivative;
	}
	
	inline void spline_cubic_smooth_end_value_example(std::ostream& out_stream = std::cout)
	{
		if (&out_stream != &std::cout) out_stream << "Data set" << ',' << "# points" << ',' << "Deviation" << ',' << "Time (ms)" << ',' << "Lambda" << std::endl;
		for (auto i = 0; i < spline_test_data_file_names.size(); ++i)
		{
			std::string name = "SplineTestData/" + spline_test_data_file_names[i];
			std::vector<tobject> x;
			std::vector<tobject> y;
			std::vector<tobject> ref;
			read_and_prepare_data(name, x, y, ref);
			int points_num = x.size();
			tobject valueLeft = ref[0];
			tobject valueRight = ref[points_num - 1];
			spline_cubic_smooth_end_value spline(valueLeft, valueRight);
			spline.load(x, y);
			auto clock_time = clock();
			spline.calculate();
			clock_time = clock() - clock_time;
			std::vector<tobject> y_f = spline.fitted_values();
			tobject dev = delta(y_f, ref);
			tobject lambda = spline.get_lambda();

			if (&out_stream != &std::cout)
			{
				std::ofstream out("output/SplineTestData/EndValue/" + spline_test_data_file_names.at(i) + "_end_value_output.csv");
				print_csv(x, y, ref, y_f, out);
				out.close();

				out_stream << spline_test_data_file_names[i] << ',' << convertToStr(&points_num) << ',' << dev << ',' << clock_time << ',' << convertToStr(&lambda) << std::endl;
			}
			std::cout << std::endl << spline_test_data_file_names[i] << std::endl
				<< "# points: " + convertToStr(&points_num) << std::endl
				<< "Deviation: " << dev << std::endl
				<< "Lambda: " << convertToStr(&lambda) << std::endl
				<< "Time: " << clock_time << " ms" << std::endl;
		}
	}
	
	inline void spline_cubic_smooth_slope_example(std::ostream& out_stream = std::cout)
	{
		if (&out_stream != &std::cout) out_stream << "Data set" << ',' << "# points" << ',' << "Deviation" << ',' << "Time (ms)" << ',' << "Lambda" << std::endl;
		for (auto i = 0; i < spline_test_data_file_names.size(); ++i)
		{
			std::string name = "SplineTestData/" + spline_test_data_file_names[i];
			std::vector<tobject> x;
			std::vector<tobject> y;
			std::vector<tobject> ref;
			read_and_prepare_data(name, x, y, ref);
			int points_num = x.size();
			tobject slopeLeft = (ref[1] - ref[0]) / (x[1] - x[0]);
			tobject slopeRight = (ref[points_num - 1] - ref[points_num - 2]) / (x[points_num - 1] - x[points_num - 2]);
			spline_cubic_smooth_slope spline(slopeLeft, slopeRight);
			spline.load(x, y);
			auto clock_time = clock();
			spline.calculate();
			clock_time = clock() - clock_time;
			std::vector<tobject> y_f = spline.fitted_values();
			tobject dev = delta(y_f, ref);
			tobject lambda = spline.get_lambda();

			if (&out_stream != &std::cout)
			{
				std::ofstream out("output/SplineTestData/Slope/" + spline_test_data_file_names.at(i) + "_slope_output.csv");
				print_csv(x, y, ref, y_f, out);
				out.close();

				out_stream << spline_test_data_file_names[i] << ',' << convertToStr(&points_num) << ',' << dev << ',' << clock_time << ',' << convertToStr(&lambda) << std::endl;
			}
			std::cout << std::endl << spline_test_data_file_names[i] << std::endl
				<< "# points: " + convertToStr(&points_num) << std::endl
				<< "Deviation: " << dev << std::endl
				<< "Lambda: " << convertToStr(&lambda) << std::endl
				<< "Time: " << clock_time << " ms" << std::endl;
		}
	}
	
	inline void spline_cubic_smooth_slope_value_example(std::ostream& out_stream = std::cout)
	{
		if (&out_stream != &std::cout) out_stream << "Data set" << ',' << "# points" << ',' << "Deviation" << ',' << "Time (ms)" << ',' << "Lambda" << std::endl;
		for (auto i = 0; i < spline_test_data_file_names.size(); ++i)
		{
			std::string name = "SplineTestData/" + spline_test_data_file_names[i];
			std::vector<tobject> x;
			std::vector<tobject> y;
			std::vector<tobject> ref;
			read_and_prepare_data(name, x, y, ref);
			int points_num = x.size();
			tobject valueLeft = ref[0];
			tobject valueRight = ref[points_num - 1];
			tobject slopeLeft = (ref[1] - ref[0]) / (x[1] - x[0]);
			tobject slopeRight = (ref[points_num - 1] - ref[points_num - 2]) / (x[points_num - 1] - x[points_num - 2]);
			spline_cubic_smooth_slope_value spline(valueLeft, valueRight, slopeLeft, slopeRight);
			spline.load(x, y);
			auto clock_time = clock();
			spline.calculate();
			clock_time = clock() - clock_time;
			std::vector<tobject> y_f = spline.fitted_values();
			tobject dev = delta(y_f, ref);
			tobject lambda = spline.get_lambda();

			if (&out_stream != &std::cout)
			{
				std::ofstream out("output/SplineTestData/SlopeValue/" + spline_test_data_file_names.at(i) + "_slope_value_output.csv");
				print_csv(x, y, ref, y_f, out);
				out.close();

				out_stream << spline_test_data_file_names[i] << ',' << convertToStr(&points_num) << ',' << dev << ',' << clock_time << ',' << convertToStr(&lambda) << std::endl;
			}
			std::cout << std::endl << spline_test_data_file_names[i] << std::endl
				<< "# points: " + convertToStr(&points_num) << std::endl
				<< "Deviation: " << dev << std::endl
				<< "Lambda: " << convertToStr(&lambda) << std::endl
				<< "Time: " << clock_time << " ms" << std::endl;
		}
	}

	inline void spline_cubic_smooth_example(std::ostream& out_stream = std::cout)
	{
		if (&out_stream != &std::cout) out_stream << "Data set" << ',' << "# points" << ',' << "Deviation" << ',' << "Time (ms)" << ',' << "Lambda" << std::endl;
		for (auto i = 0; i < spline_test_data_file_names.size(); ++i)
		{
			std::string name = "SplineTestData/" + spline_test_data_file_names[i];
			std::vector<tobject> x;
			std::vector<tobject> y;
			std::vector<tobject> ref;
			read_and_prepare_data(name, x, y, ref);
			int points_num = x.size();
			//tobject valueLeft = ref[0];
			//tobject valueRight = ref[points_num - 1];
			//tobject slopeLeft = (ref[1] - ref[0]) / (x[1] - x[0]);
			//tobject slopeRight = (ref[points_num - 1] - ref[points_num - 2]) / (x[points_num - 1] - x[points_num - 2]);
			smooth_cubic_spline spline;
			spline.load(x, y);
			auto clock_time = clock();
			spline.calculate();
			clock_time = clock() - clock_time;
			std::vector<tobject> y_f = spline.fitted_values();
			tobject dev = delta(y_f, ref);
			tobject lambda = spline.get_lambda();

			if (&out_stream != &std::cout)
			{
				std::ofstream out("output/SplineTestData/" + spline_test_data_file_names.at(i) + "_output.csv");
				print_csv(x, y, ref, y_f, out);
				out.close();

				out_stream << spline_test_data_file_names[i] << ',' << convertToStr(&points_num) << ',' << dev << ',' << clock_time << ',' << convertToStr(&lambda) << std::endl;
			}
			std::cout << std::endl << spline_test_data_file_names[i] << std::endl
				<< "# points: " + convertToStr(&points_num) << std::endl
				<< "Deviation: " << dev << std::endl
				<< "Lambda: " << convertToStr(&lambda) << std::endl
				<< "Time: " << clock_time << " ms" << std::endl;
		}
	}

    inline void smooth_cubic_spline_example(std::ostream& out_stream = std::cout)
    {
        const int n = 10;
		std::vector<tobject> x(n + 1);
        std::vector<tobject> y(n + 1);
        boost::mt19937 rng;
        boost::normal_distribution<double> nd(0.0, 0);
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<>> var_nor(rng, nd);
        for (int i = 0; i < n + 1; i++)
        {
            x[i] = (double)i / n;
            tobject noise = var_nor();
            y[i] = func(x[i]) + noise;
        }
        out_stream << "x " << x << "\n";
        out_stream << "y " << y << "\n";
        smooth_cubic_spline spline;
        spline.load(x, y);
        spline.calculate();
        std::vector<tobject> y_f = spline.fitted_values();
        out_stream << "y_fitted " << y_f << "\n";
        tobject dev = deviation(&spline, y_f, func);
        out_stream << std::endl << "deviation between estimation and reference: " << dev << std::endl << std::endl;
    }

	// Auxiliary method for generating noisy data.
	inline void data(int n, std::vector<double>& x, std::vector<tobject>& y, std::function<tobject(double)> function, boost::normal_distribution<double> nd)
	{
		x.resize(n + 1);
		y.resize(n + 1);
		boost::mt19937 rng;
		boost::variate_generator<boost::mt19937&, boost::normal_distribution<>> var_nor(rng, nd);
		for (int i = 0; i < n + 1; i++)
		{
			x[i] = (double)i / n;
			tobject noise = var_nor();
			y[i] = function(x[i]) + noise;
		}
	}

   // // Test for differentiation of smooth cubic spline.
   // inline void smooth_cubic_spline_example_with_tobject(std::ostream& out_stream = std::cout)
   // {
   //     const int n = 9;
   //     double step = 1e-3;
   //     std::vector<double> x;
   //     std::vector<tobject> y;
   //     boost::normal_distribution<double> nd(0.0, 5);
   //     data(n, x, y, func, nd);
   //     smooth_cubic_spline spline;
   //     tcontext<tobject> spline_context;
   //     for (int i = 0; i < y.size(); i++)
   //     {
   //         std::stringstream ss;
   //         ss.operator<<(i);
   //         std::string str("y_noised_");
   //         ss >> str;
   //         spline_context.independent(y[i], str);
   //     }
   //     std::vector<tobject> y_context(n + 1);
   //     for (int i = 0; i < n + 1; i++)
   //         y_context[i] = spline_context[i];
   //     spline.load(x, y_context);
   //     spline.calculate();
   //     std::vector<tobject> y_f = spline.fitted_values();
   //     for (int i = 0; i < y.size(); i++)
   //     {
   //         std::stringstream ss;
			//ss.operator<<(i);
   //         std::string str("y_fitted_");
   //         ss >> str;
   //         spline_context.dependent(y_f[i], str);
   //     }
   //     std::vector<double> forward;
   //     std::vector<double> reverse;
   //     std::vector<double> dx(y.size());
   //     long finite_difference_time = -std::clock();
   //     std::vector<std::vector<tobject>> finite_difference = finite_difference_derivatives(x, y, step);
   //     finite_difference_time += std::clock();
   //     long forward_time = 0;
   //     long reverse_time = 0;
   //     out_stream << "Differentiating fitted y values by input y values... \n\n";
   //     for (int i = 0; i < y.size(); i++)
   //     {
   //         for (int j = 0; j < y.size(); j++)
   //             dx[j] = 0;
   //         dx[i] = 1;
   //         forward_time -= std::clock();
   //         forward = spline_context.forward(1, dx);
   //         forward_time += std::clock();
   //         reverse_time -= std::clock();
   //         reverse = spline_context.reverse(1, dx);
   //         reverse_time += std::clock();
   //         out_stream << "Forward derivatives = " << forward << "..." << std::endl;
   //         out_stream << "Reverse derivatives = " << reverse << "..." << std::endl;
   //         out_stream << "Finite-difference derivatives = " << finite_difference[i] << "..." << std::endl;
   //     }
   //     out_stream << "Forward calculation took (ms): " << forward_time / (double)(CLOCKS_PER_SEC)* 1000 << std::endl;
   //     out_stream << "Reverse calculation took (ms): " << reverse_time / (double)(CLOCKS_PER_SEC)* 1000 << std::endl;
   //     out_stream << "Finite-difference calculation took (ms): " << finite_difference_time / (double)(CLOCKS_PER_SEC)* 1000 << std::endl;
   // }

	// Test for differentiation of smooth cubic spline.
	inline void spline_cubic_smooth_slope_value_diff_example(std::ostream& out_stream = std::cout)
	{
		const int point_number = 5;
		double step = 1e-3;
		const bool print_performance = true;
		// Array with file names.
		std::vector<std::string> names =
		{
			"1_RiskReversal",
			"2_Straddle",
			"3_Strangle",
			"AMCDataBear1D1000",
			"AMCDataBull1D1000",
			"AMCDataButterfly1D1000",
			"DownAndInCallBarrierTest1000",
			"UpAndInCallBarrierTest1000",
			"UpAndOutCallBarrierTest1000"
		};
		int files_number = names.size();
		out_stream << "\nTest of Adjoint_Y_slope_value_example \n";
		if (&out_stream != &std::cout) std::cout << "\nTest of Adjoint_Y_slope_value_example \n";
		for (int nn = 0; nn < files_number; nn++)
		{
			out_stream << "\n Test for " << names[nn] << "\n";
			if (&out_stream != &std::cout) std::cout << "\n Test for " << names[nn] << "\n";
			std::vector<tobject> x;
			std::vector<tobject> y;
			std::vector<tobject> ref;
			std::string name = "SplineTestData/" + names[nn];
			read_and_prepare_data(name, x, y, ref);

			int n = y.size();
			int n_mult = n / (point_number + 1); // if end or slope value here we can't change 0 and n - 1 points
			std::vector<tobject> X = { cl::tapescript::pack_vec(y), cl::tapescript::pack_vec(x) };

			tape_start(X);

			tobject valueLeft = ref[0];
			tobject valueRight = ref[n - 1];
			tobject slopeLeft = (ref[1] - ref[0]) / (x[1] - x[0]);
			tobject slopeRight = (ref[n - 1] - ref[n - 2]) / (x[n - 1] - x[n - 2]);
			spline_cubic_smooth_slope_value spline(valueLeft, valueRight, slopeLeft, slopeRight);
			spline.load(cl::tapescript::unpack_vec(X[1], x.size()), cl::tapescript::unpack_vec(X[0], y.size()));
			spline.set_lambda(estimate_rigidity(y)); // cl::tapescript::unpack_vec(X[0], y.size()))); // if want lambda Independent
			spline.calculate();
			tobject estim = cl::tapescript::pack_vec(spline.fitted_values());

			std::vector<tobject> Y = { estim };

			tfunc<tvalue> f(X, Y);

			// Finite difference derivatives.
			long finite_difference_time = -std::clock();
			std::vector<std::vector<tobject>> derivative = finite_difference_derivatives(spline, x, y, point_number, step);
			finite_difference_time += std::clock();
			//if (print_performance)
			//	out_stream << " Calculation time of Finite difference (ms): " << finite_difference_time / point_number << "\n\n";

			std::vector<cl::tvalue> forward;
			std::vector<cl::tvalue> reverse;
			cl::tvalue zero;
			zero.resize(n);

			long f_time = 0;
			long r_time = 0;
			long forward_time = 0;
			long reverse_time = 0;
			out_stream << "\n  Differentiating fitted y values with respect to input " << (diff_with_respect_to_y ? "y" : "x") << " values... \n\n";
			std::vector<std::vector<tobject>> forward_derivatives(point_number);
			std::vector<std::vector<tobject>> reverse_derivatives(point_number);
			for (int i = 0; i < point_number; i++)
			{
				forward_derivatives[i].resize(point_number);
				reverse_derivatives[i].resize(point_number);
			}
			for (int i = 0; i < point_number; i++)
			{
				std::vector<cl::tvalue> dy = { zero, zero };
				dy[diff_with_respect_to_y ? 0 : 1][(i + 1) * n_mult] = 1;
				forward_derivatives[i].resize(point_number);
				f_time -= std::clock();
				forward = f.forward(1, dy);
				f_time += std::clock();
				forward_time += f_time;
				std::vector<cl::tvalue> dy0(1, zero);
				dy0[0][(i + 1) * n_mult] = 1;
				r_time -= std::clock();
				reverse = f.reverse(1, dy0);
				r_time += std::clock();
				reverse_time += r_time;
				for (int j = 0; j < point_number; j++)
				{
					forward_derivatives[i][j] = forward[0][(j + 1) * n_mult];
					reverse_derivatives[j][i] = reverse[diff_with_respect_to_y ? 0 : 1][(j + 1) * n_mult];
				}
			}
			std::cout << std::endl;
			out_stream << "\tPoints number = " << point_number << std::endl;
			out_stream << "\tLambda = " << spline.get_lambda() << std::endl << std::endl;
			for (int i = 0; i < point_number; i++)
			{
				out_stream << "\tAt point " << (i + 1) * n_mult + 1<< " of " << n << std::endl << std::endl;
				out_stream << "\tForward derivatives = " << forward_derivatives[i] << std::endl;
				out_stream << "\tReverse derivatives = " << reverse_derivatives[i] << std::endl;
				out_stream << "\tFinite-difference derivatives = " << derivative[i] << std::endl << std::endl;
				for (int j = 0; j < point_number; j++)
				{
					//out_stream << "\t\t|Forward - Reverse| = " << std::abs(forward_derivatives[i][j] - reverse_derivatives[i][j]) << std::endl;
					out_stream << "\t\t|Forward - Finite-difference| = " << (std::abs(forward_derivatives[i][j] - derivative[i][j])) << std::endl;
					out_stream << "\t\t|Reverse - Finite-difference| = " << (std::abs(reverse_derivatives[i][j] - derivative[i][j])) << std::endl << std::endl;
				}
			}
			if (print_performance)
			{
				out_stream << " Calculation time of Forward (ms): " << forward_time / point_number << "\n";
				out_stream << " Calculation time of Reverse (ms): " << reverse_time / point_number << "\n";
				out_stream << " Calculation time of Finite difference (ms): " << finite_difference_time / point_number << "\n\n";
			}
			out_stream << " Test for " << names[nn] << " completed.\n";
		}
		out_stream << "Test completed. \n\n";
		if (&out_stream != &std::cout) std::cout << "Test completed. \n\n";
	}

    inline void cubic_spline_examples()
    {
        std::ofstream of("output/smooth_cubic_spline_example_output.txt");
        cl::tape_serializer<cl::tvalue> serializer(of);
        serializer.precision(3);

        smooth_cubic_spline_example(serializer);
        //smooth_cubic_spline_example_with_tdouble(serializer);
    }
}
#endif // cl_tape_examples_impl_smooth_cubic_spline_example_hpp