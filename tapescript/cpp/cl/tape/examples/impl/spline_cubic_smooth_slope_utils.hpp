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

#ifndef cl_tape_examples_impl_spline_cubic_smooth_slope_utils_hpp
#define cl_tape_examples_impl_spline_cubic_smooth_slope_utils_hpp

#include <cl/tape/tape.hpp>
#include "impl/utils.hpp"
#include "impl/smooth_cubic_spline_utils.hpp"

namespace cl
{
	/// <summary>Math utility class for spline_cubic_smooth_slope.</summary>
	class spline_cubic_smooth_slope_utils : public smooth_cubic_spline_utils
	{
	private:
		tobject slopeLeft_;
		tobject slopeRight_;
		std::vector<tobject> weight;
	protected:
		///<summary> Initialise right side of equation. </summary>
		void calculate_q() override
		{
			auxiliaries[4].resize(data_number);
			auxiliaries[4][0] = 3 * (y[1] - y[0]) / auxiliaries[1][0] - 3 * slopeLeft_;
			for (int i = 1; i < data_number - 1; ++i)
				auxiliaries[4][i] = 3 * (y[i + 1] - y[i]) / auxiliaries[1][i] - 3 * (y[i] - y[i - 1]) / auxiliaries[1][i - 1];
			auxiliaries[4][data_number - 1] = -3 * (y[data_number - 1] - y[data_number - 2]) / auxiliaries[1][data_number - 2] + 3 * slopeRight_;
		}

		///<summary> Prepare system of linear equations: calculate all non-zero matrix elements.</summary>
		void prepare_linear_system(tobject mu) override
		{
			int n = data_number;
			std::vector<std::vector<tobject>>& aux = auxiliaries;
			//aux[0].resize(n);
			//aux[3].resize(n);
			uvw_.resize(3);
			std::vector<tobject>& u = uvw_[0];
			std::vector<tobject>& v = uvw_[1];
			std::vector<tobject>& w = uvw_[2];
			// Diagonal, first and second upper diagonal vectors u, v, w, respectively.
			u.resize(n);
			v.resize(n - 1);
			w.resize(n - 2);
			u[0] = mu * (aux[3][0] * aux[3][0] * weight[0] + aux[2][0] * aux[2][0] * weight[1]) + aux[0][0];
			v[0] = mu * aux[2][0] * (aux[3][0] * weight[0] + aux[3][1] * weight[1]) + aux[1][0];
			w[0] = mu*aux[2][0] * aux[2][1] * weight[1];
			for (int i = 1; i < n - 2; i++)
			{
				u[i] = mu * (aux[2][i - 1] * aux[2][i - 1] * weight[i - 1] + aux[3][i] * aux[3][i] * weight[i]
							 + aux[2][i] * aux[2][i] * weight[i + 1]) + aux[0][i];
				v[i] = mu * aux[2][i] * (aux[3][i] * weight[i] + aux[3][i + 1] * weight[i + 1]) + aux[1][i];
				w[i] = mu * aux[2][i] * aux[2][i + 1] * weight[i + 1];
			}
			u[n - 2] = mu * (aux[2][n - 3] * aux[2][n - 3] * weight[n - 3]
										 + aux[3][n - 2] * aux[3][n - 2] * weight[n - 2]
										 + aux[2][n - 2] * aux[2][n - 2] * weight[n - 1]) + aux[0][n - 2];
			u[n - 1] = mu * (aux[2][n - 2] * aux[2][n - 2] * weight[n - 2]
										 + aux[3][n - 1] * aux[3][n - 1] * weight[n - 1]) + aux[0][n - 1];
			v[n - 2] = mu * aux[2][n - 2] * (aux[3][n - 2] * weight[n - 2]
																 + aux[3][n - 1] * weight[n - 1]) + aux[1][n - 2];
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
			D[0] = y[0] - mu * weight[0] * (aux[3][0] * aux[4][0] + aux[2][0] * aux[4][1]);
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
			D[n - 1] = y[n - 1] - mu * weight[n - 1] * (aux[2][n - 2] * aux[4][n - 2] + aux[3][n - 1] * aux[4][n - 1]);
			coeff[0][n - 1] = 0;
			coeff[2][n - 1] = slopeRight_;
		}

		void calculate_helper_variables(const std::vector<tobject>& x0, std::vector<std::vector<tobject>>& aux) override
		{
			aux.resize(5);
			int n = x0.size();
			aux[0].resize(n);
			aux[1].resize(n - 1);
			aux[2].resize(n - 1);
			aux[3].resize(n);
			aux[1][0] = x0[1] - x0[0];
			aux[2][0] = 3 / aux[1][0];
			aux[0][0] = 2 * aux[1][0];
			aux[3][0] = -aux[2][0];
			for (int i = 1; i < n - 1; ++i)
			{
				aux[1][i] = x0[i + 1] - x0[i];
				aux[2][i] = 3 / aux[1][i];
				aux[0][i] = 2 * (aux[1][i - 1] + aux[1][i]);
				aux[3][i] = -(aux[2][i - 1] + aux[2][i]);
			}
			aux[0][n - 1] = 2 * aux[1][n - 2];
			aux[3][n - 1] = -aux[2][n - 2];
		}
	public:
		spline_cubic_smooth_slope_utils(tobject slopeLeft, tobject slopeRight) : slopeLeft_(slopeLeft), slopeRight_(slopeRight) {}
		/// <summary> Spline initialisation from X, Y input arrays and boundary conditions (slopes at end points)
		/// X must be in ascending order. </summary>
		void load(const std::vector<tobject>& x0, const std::vector<tobject>& y0, const std::vector<tobject>& weight0, tobject slopeLeft, tobject slopeRight)
		{
			weight.clear();
			weight = weight0;

			load(x0, y0, slopeLeft, slopeRight);
		}
		void load(const std::vector<tobject>& x0, const std::vector<tobject>& y0, tobject slopeLeft, tobject slopeRight)
		{
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
#endif // cl_tape_examples_impl_spline_cubic_smooth_slope_utils_hpp