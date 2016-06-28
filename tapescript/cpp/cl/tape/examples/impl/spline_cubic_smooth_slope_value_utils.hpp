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

#ifndef cl_tape_examples_impl_spline_cubic_smooth_slope_value_utils_hpp
#define cl_tape_examples_impl_spline_cubic_smooth_slope_value_utils_hpp

#include <cl/tape/tape.hpp>
#include "impl/utils.hpp"
#include "impl/smooth_cubic_spline_utils.hpp"
#include "impl/spline_cubic_smooth_end_value_utils.hpp"

namespace cl
{
	/// <summary>Math utility class for spline_cubic_smooth_slope_value.</summary>
	class spline_cubic_smooth_slope_value_utils : public spline_cubic_smooth_end_value_utils
	{
	private:
		tobject slopeLeft_;
		tobject slopeRight_;
	protected:
		///<summary> Initialise right side of equation. </summary>
		void calculate_q() override
		{
			spline_cubic_smooth_end_value_utils::calculate_q();
			auxiliaries[4][0] += -3 * slopeLeft_;
			auxiliaries[4][data_number - 1] += 3 * slopeRight_;
		}
		
		///<summary> Calculate all spline coefficients.</summary>
		void find_all_coefficients(tobject mu, std::vector<tobject>& D) override
		{
			int n = data_number;
			D.resize(n); // n - 1 -> n
			for (int i = 0; i < 3; i++)
				coefficients[i].resize(n);
			std::vector<std::vector<tobject>>& aux = auxiliaries;
			std::vector<std::vector<tobject>>& coeff = coefficients;
			coeff[1][0] = aux[4][0];
			D[0] = valueLeft_;
			coeff[0][0] = (aux[4][1] - aux[4][0]) / (3 * aux[1][0]);
			coeff[2][0] = slopeLeft_;
			for (int i = 1; i < n - 1; ++i)
			{
				coeff[1][i] = aux[4][i];
				D[i] = y[i] - mu * weight[i] * (aux[2][i - 1] * aux[4][i - 1] + aux[3][i] * aux[4][i] + aux[2][i] * aux[4][i + 1]);
				coeff[0][i] = (aux[4][i + 1] - aux[4][i]) / (3 * aux[1][i]);
				coeff[2][i] = (aux[4][i] + aux[4][i - 1]) * aux[1][i - 1] + coeff[2][i - 1];
			}
			coeff[1][n - 1] = aux[4][n - 1];
			D[n - 1] = valueRight_;
			coeff[0][n - 1] = 0;
			coeff[2][n - 1] = slopeRight_;
		}
	public:
		spline_cubic_smooth_slope_value_utils(tobject valueLeft, tobject valueRight, tobject slopeLeft, tobject slopeRight) : slopeLeft_(slopeLeft), slopeRight_(slopeRight)
		{
			valueLeft_ = valueLeft;
			valueRight_ = valueRight;
		}
		/// <summary> Spline initialisation from X, Y input arrays and boundary conditions (values at end points)
		/// X must be in ascending order. </summary>
		void load(const std::vector<tobject>& x0, const std::vector<tobject>& y0, const std::vector<tobject>& weight0, tobject valueLeft, tobject valueRight, tobject slopeLeft, tobject slopeRight)
		{
			weight.clear();
			weight = weight0;

			load(x0, y0, valueLeft, valueRight, slopeLeft, slopeRight);
		}
		void load(const std::vector<tobject>& x0, const std::vector<tobject>& y0, tobject valueLeft, tobject valueRight, tobject slopeLeft, tobject slopeRight)
		{
			valueLeft_ = valueLeft;
			valueRight_ = valueRight;
			slopeLeft_ = slopeLeft;
			slopeRight_ = slopeRight;

			load(x0, y0);
		}
		void load(const std::vector<tobject>& x0, const std::vector<tobject>& y0) override
		{
			// if weight data size don't have the same items as X
			if (weight.size() < x0.size())
			for (int i = weight.size(); i < x0.size(); ++i)
				weight.emplace_back(1);
			if (weight.size() > x0.size())
			for (int i = weight.size(); i > x0.size(); --i)
				weight.pop_back();

			smooth_cubic_spline_utils::load(x0, y0);
		}
	};
}
#endif // cl_tape_examples_impl_spline_cubic_smooth_slope_value_utils_hpp