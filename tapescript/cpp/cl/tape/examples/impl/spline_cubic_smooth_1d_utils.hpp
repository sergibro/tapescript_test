/*
Copyright (C) 2013-present CompatibL. All rights reserved.

This file is part of ModVal Engine (the "Software"), a model validation
library available from:

http://git.modval.org (source)
http://www.modval.org (documentation)

The Software is distributed under multiple licenses. This distribution
is under the terms of the ModVal.org license (the "License").
You may obtain a copy of the License at:

http://www.modval.org/about/license/

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF TITLE,
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.

The Software is licensed for regulatory and internal model validation
use only, subject to the License and the inclusion of this copyright
notice. The use of the Software or its derivative works, in whole or
in part, in trading, risk management, consulting, or commercial software,
or its redistribution in modified or unmodified form is prohibited
without prior written permission.
*/

#ifndef cl_tape_examples_impl_spline_cubic_smooth_1d_utils_hpp
#define cl_tape_examples_impl_spline_cubic_smooth_1d_utils_hpp

#include <cl/tape/tape.hpp>
#include "impl/utils.hpp"

namespace cl
{
	/// <summary>Math utility class for ClSplineBase.</summary>
	class spline_cubic_smooth_1d_utils
	{
	protected:
		// Helper variables needed for fast and compact calculations.
		// (their one-letter names exactly follow Pollock's article).
		std::vector<tdouble> p_;
		std::vector<tdouble> h_;
		std::vector<tdouble> r_;
		std::vector<tdouble> f_;
		std::vector<tdouble> q_;

		// Dimension of the system of linear equations (number of unknown b_i).
		int equation_size_;

		///<summary> Initialise spline helper variables. </summary>
		void calculate_helper_variables()
		{
			p_.resize(points_number_);
			h_.resize(points_number_ - 1);
			r_.resize(points_number_ - 1);
			f_.resize(points_number_);
			h_[0] = spline_sections_[1].X - spline_sections_[0].X;
			r_[0] = 3 / h_[0];
			p_[0] = 2 * h_[0];
			f_[0] = -r_[0];
			for (int i = 1; i < points_number_ - 1; ++i)
			{
				h_[i] = spline_sections_[i + 1].X - spline_sections_[i].X;
				r_[i] = 3 / h_[i];
				p_[i] = 2 * (h_[i - 1] + h_[i]);
				f_[i] = -(r_[i - 1] + r_[i]);
			}
			p_[points_number_ - 1] = 2 * h_[points_number_ - 2];
			f_[points_number_ - 1] = -r_[points_number_ - 2];
			calculate_Q();
		}

		///<summary> Initialise right side of equation. </summary>
		virtual void calculate_Q()
		{
			q_.resize(points_number_ - 2);
			for (int i = 1; i < points_number_ - 1; ++i)
				q_[i - 1] = 3 * (spline_sections_[i + 1].Y - spline_sections_[i].Y) / h_[i] - 3 * (spline_sections_[i].Y - spline_sections_[i - 1].Y) / h_[i - 1];
		}

		///<summary> Prepare system of linear equations: calculate all non-zero matrix elements.</summary>
		virtual void prepare_linear_system(tdouble mu, std::vector<tdouble> &u, std::vector<tdouble> &v, std::vector<tdouble> &w)
		{
			// Diagonal, first and second upper diagonal vectors u, v, w, respectively.
			u.resize(equation_size_);
			v.resize(equation_size_ - 1);
			w.resize(equation_size_ - 2);
			for (int i = 0; i < equation_size_ - 2; i++)
			{
				u[i] = mu * (r_[i] * r_[i] * spline_sections_[i].Weight + f_[i + 1] * f_[i + 1] * spline_sections_[i + 1].Weight
							 + r_[i + 1] * r_[i + 1] * spline_sections_[i + 2].Weight) + p_[i + 1];
				v[i] = mu * (r_[i + 1] * (f_[i + 1] * spline_sections_[i + 1].Weight + f_[i + 2] * spline_sections_[i + 2].Weight)) + h_[i + 1];
				w[i] = mu * r_[i + 1] * r_[i + 2] * spline_sections_[i + 2].Weight;
			}
			u[equation_size_ - 2] = mu * (r_[equation_size_ - 2] * r_[equation_size_ - 2] * spline_sections_[equation_size_ - 2].Weight +
										 f_[equation_size_ - 1] * f_[equation_size_ - 1] * spline_sections_[equation_size_ - 1].Weight +
										 r_[equation_size_ - 1] * r_[equation_size_ - 1] * spline_sections_[equation_size_].Weight) + p_[equation_size_ - 1];
			u[equation_size_ - 1] = mu * (r_[equation_size_ - 1] * r_[equation_size_ - 1] * spline_sections_[equation_size_ - 1].Weight +
										 f_[equation_size_] * f_[equation_size_] * spline_sections_[equation_size_].Weight +
										 r_[equation_size_] * r_[equation_size_] * spline_sections_[equation_size_ + 1].Weight) + p_[equation_size_];
			v[equation_size_ - 2] = mu * r_[equation_size_ - 1] * (f_[equation_size_ - 1] * spline_sections_[equation_size_ - 1].Weight +
																 f_[equation_size_] * spline_sections_[equation_size_].Weight) + h_[equation_size_ - 1];
		}

		///<summary> Calculate spline B coefficients 
		/// (solution of linear system of equations).</summary>
		void find_B_coefficients(std::vector<tdouble> &u, std::vector<tdouble> &v, std::vector<tdouble> &w)
		{
			// Factorization procedure.
			w[0] /= u[0];
			v[0] /= u[0];
			u[1] -= u[0] * v[0] * v[0];
			v[1] = (v[1] - u[0] * v[0] * w[0]) / u[1];
			w[1] /= u[1];
			for (int i = 3; i <= equation_size_ - 2; i++)
			{
				u[i - 1] -= (u[i - 3] * w[i - 3] * w[i - 3] + u[i - 2] * v[i - 2] * v[i - 2]);
				v[i - 1] = (v[i - 1] - u[i - 2] * v[i - 2] * w[i - 2]) / u[i - 1];
				w[i - 1] /= u[i - 1];
			}
			u[equation_size_ - 2] -= (u[equation_size_ - 4] * w[equation_size_ - 4] * w[equation_size_ - 4] +
									 u[equation_size_ - 3] * v[equation_size_ - 3] * v[equation_size_ - 3]);
			v[equation_size_ - 2] = (v[equation_size_ - 2] - u[equation_size_ - 3] * v[equation_size_ - 3] * w[equation_size_ - 3]) / u[equation_size_ - 2];
			u[equation_size_ - 1] -= (u[equation_size_ - 3] * w[equation_size_ - 3] * w[equation_size_ - 3] +
									 u[equation_size_ - 2] * v[equation_size_ - 2] * v[equation_size_ - 2]);

			// Forward substitution.
			q_[1] -= v[0] * q_[0];
			for (int i = 2; i < equation_size_; ++i)
				q_[i] -= (v[i - 1] * q_[i - 1] + w[i - 2] * q_[i - 2]);
			for (int i = 0; i < equation_size_; ++i)
				q_[i] /= u[i];

			// Backward substitution.
			q_[equation_size_ - 2] = q_[equation_size_ - 2] - v[equation_size_ - 2] * q_[equation_size_ - 1];
			for (int i = equation_size_ - 3; i >= 0; --i)
				q_[i] -= (v[i] * q_[i + 1] + w[i] * q_[i + 2]);
		}

		///<summary> Calculate all spline coefficients.</summary>
		virtual void find_all_coefficients(tdouble mu)
		{
			spline_sections_[0].D = spline_sections_[0].Y - mu * r_[0] * q_[0] * spline_sections_[0].Weight;
			spline_sections_[1].D = spline_sections_[1].Y - mu * spline_sections_[0].Weight * (f_[1] * q_[0] + r_[1] * q_[1]);
			spline_sections_[0].A = q_[0] / (3 * h_[0]);
			spline_sections_[0].B = 0.0;
			spline_sections_[0].C = (spline_sections_[1].D - spline_sections_[0].D) / h_[0] - q_[0] * h_[0] / 3;
			spline_sections_[1].A = (q_[1] - q_[0]) / (3 * h_[1]);
			spline_sections_[1].B = q_[0];
			spline_sections_[1].C = q_[0] * h_[0] + spline_sections_[0].C;
			for (int i = 2; i < points_number_ - 2; ++i)
			{
				spline_sections_[i].A = (q_[i] - q_[i - 1]) / (3 * h_[i]);
				spline_sections_[i].B = q_[i - 1];
				spline_sections_[i].C = (q_[i - 1] + q_[i - 2]) * h_[i - 1] + spline_sections_[i - 1].C;
				spline_sections_[i].D = spline_sections_[i].Y - mu * spline_sections_[i].Weight * (r_[i - 1] * q_[i - 2] + f_[i] * q_[i - 1] + r_[i] * q_[i]);
			}
			spline_sections_[points_number_ - 2].A = (0 - q_[points_number_ - 3]) / (3 * h_[points_number_ - 2]);
			spline_sections_[points_number_ - 2].B = q_[points_number_ - 3];
			spline_sections_[points_number_ - 2].C = (q_[points_number_ - 3] + q_[points_number_ - 4]) * h_[points_number_ - 3] + spline_sections_[points_number_ - 3].C;
			spline_sections_[points_number_ - 2].D = spline_sections_[points_number_ - 2].Y - mu * spline_sections_[points_number_ - 2].Weight
				* (r_[points_number_ - 3] * q_[points_number_ - 4] + f_[points_number_ - 2] * q_[points_number_ - 3]);
		}
	public:
		// Number of data point and spline sections.
		int points_number_;

		// Array of spline sections.
		std::vector<spline_cubic_selection> spline_sections_;

		/// <summary> Spline initialisation from X, Y input arrays.
		/// u is array of uncertainties.
		/// X must be in ascending order. </summary>
		virtual void load(std::vector<weighted_point> data, spline_1d_builder  parameters) //parameters = NULL
		{
			// Number of points.
			points_number_ = data.size();

			// Number of unknown b_i.
			equation_size_ = points_number_ - 2;

			// Initialise spline sections (input points are (x_0, y_0), ..., (x_n, y_n)).
			spline_sections_.resize(points_number_);
			for (int i = 0; i < points_number_; ++i)
			{
				spline_sections_[i].reload(data[i].x, data[i].y, 1); //X = data[i].X[0], Y = data[i].Value, Weight = 1
				// TODO Fix NaN problem with weigths.
				//spline_sections_[i].Weight = data[i].Weight;
			}

			// Calculate spline helper variables.
			calculate_helper_variables();
		}

		/// <summary>Returns knot section coordinate.</summary>
		tdouble get_knot_value(int section)
		{
			return spline_sections_[section].X;
		}

		/// <summary>Returns data value in section.</summary>
		tdouble get_data_value(int section)
		{
			return spline_sections_[section].Y;
		}

		/// <summary>Get spline derivative at point x with possible extrapolation.</summary>
		tdouble value_at(tdouble x, int derivative, bool extrapolate)
		{
			// If spline can be extrapolated, check if requested x is outside its sections.
			if (extrapolate)
			{
				if (x < spline_sections_[0].X)
					return get_spline_value_in_section(0, x - spline_sections_[0].X, derivative);
				else if (x > spline_sections_[spline_sections_.size() - 1].X)
					return get_spline_value_in_section(spline_sections_.size() - 2, x - spline_sections_[spline_sections_.size() - 2].X, derivative);
			}

			// Check for requested x within spline.
			// TODO Use binary search here.
			for (int j = 0; j < spline_sections_.size() - 1; ++j)
			{
				if (x >= spline_sections_[j].X && x <= spline_sections_[j + 1].X)
					return get_spline_value_in_section(j, x - spline_sections_[j].X, derivative);
			}

			// If spline section for requested x was not found, return empty.
			throw "Value not found.";// ClEx("Value not found.");
		}

		///<summary> Return spline value in the specified section (j) at the specified internal point (x).
		/// Optionaly specify derivative (0 by default, i.e. spline value).</summary>
		tdouble get_spline_value_in_section(int j, tdouble x, int derivative)
		{
			return spline_sections_[j].value_at(x, derivative);
		}

		///// <summary>  Return ClCsvMatrix with all spline coefficients.  </summary>
		//public ClCsvMatrix GetCoefficients()
		//{
		//	ClCsvMatrix matrix = new ClCsvMatrix(spline_sections_.size() - 1, 6);
		//	for (int i = 0; i < spline_sections_.size() - 1; i++)
		//	{
		//		matrix[i, 0] = new tdouble(spline_sections_[i].X).ToVariant();
		//		matrix[i, 1] = new tdouble(spline_sections_[i + 1].X).ToVariant();
		//		matrix[i, 2] = new tdouble(spline_sections_[i].A).ToVariant();
		//		matrix[i, 3] = new tdouble(spline_sections_[i].B).ToVariant();
		//		matrix[i, 4] = new tdouble(spline_sections_[i].C).ToVariant();
		//		matrix[i, 5] = new tdouble(spline_sections_[i].D).ToVariant();
		//	}
		//	return matrix;
		//}

		///<summary> Calculate spline coefficients.</summary>
		void solve(tdouble lambda)
		{
			// Bound smoothing parameter value to avoid problems with NaN.
			const tdouble lambdaMin = 1e-6;
			const tdouble lambdaMax = 1 - 1e-6;
			if (lambda < lambdaMin)
				lambda = lambdaMin;
			else if (lambda > lambdaMax)
				lambda = lambdaMax;

			// Internal smoothing parameter: transform lambda -> mu = 2 * (1 - lambda) / (3 * lambda).
			tdouble mu = 2 * (1 - lambda) / (3 * lambda);

			// Prepare system of linear equations: calculate diagonal, first and second upper diagonal vectors u, v, w, respectively.
			std::vector<tdouble> u, v, w;
			prepare_linear_system(mu, u, v, w);

			// Find b-coefficients in splines in q_ array.
			find_B_coefficients(u, v, w);

			// Find the rest of spline parameters.
			find_all_coefficients(mu);
		}

		///<summary>Apply boundary conditions (add extra sectios) 
		/// to match specified slopes slope1 and slope2.</summary>
		void apply_extra_slope(tdouble slope1, tdouble slope2)
		{
			// Create list of existing spline sections.
			std::vector<spline_cubic_selection> list_section(spline_sections_);
			
			// Add extra sections at the beginning.
			tdouble extraStep = h_.at(0);
			std::vector<spline_cubic_selection> extra_cubic_smooth_section_begin = list_section[0].extend(slope1, extraStep);
			for (int i = 2; i >= 0; i--)
				list_section.insert(list_section.begin(), extra_cubic_smooth_section_begin.at(i));
			
			// Add extra sections at the end (this is less trivial).
			// First calculate negative extra step.
			extraStep = -h_.at(h_.size() - 1);
			spline_cubic_selection cubic_smooth_section_last = list_section[list_section.size() - 2];
			// Re-reference last spline section.
			cubic_smooth_section_last.change_start_point(list_section.at(list_section.size() - 1).X, list_section.at(list_section.size() - 1).Y);
			// Obtain needed extra sections.
			std::vector<spline_cubic_selection> extra_cubic_smooth_section_end = cubic_smooth_section_last.extend(slope2, extraStep);
			// And re-reference them all back.
			spline_cubic_selection cubicSmoothSectionTail(extra_cubic_smooth_section_end[0].X, extra_cubic_smooth_section_end[0].Y);
			extra_cubic_smooth_section_end[0].change_start_point(extra_cubic_smooth_section_end[1].X, extra_cubic_smooth_section_end[1].Y);
			extra_cubic_smooth_section_end[1].change_start_point(extra_cubic_smooth_section_end[2].X, extra_cubic_smooth_section_end[2].Y);
			extra_cubic_smooth_section_end[2].change_start_point(list_section[list_section.size() - 1].X, list_section[list_section.size() - 1].Y);
			// Finally, remove arteficial last section and add new ones.
			list_section.pop_back();
			for (int i = 2; i >= 0; i--)
				list_section.push_back(extra_cubic_smooth_section_end[i]);
			list_section.push_back(cubicSmoothSectionTail);

			// Update field spline_sections_.
			spline_sections_ = list_section;
		}

		///<summary> Bound spline knots.</summary>
		bool bound_knots(tdouble min, tdouble max, tdouble overCorrection)
		{
			// Set default argument value.
			//if (overCorrection.IsEmpty)
			//	overCorrection = 0;

			// Check knots, correct data if needed.
			bool flagModified = false;
			for (int i = 0; i < points_number_; i++)
			{
				tdouble fitValue = (i == points_number_ - 1) ? get_spline_value_in_section(i - 1, h_.at(h_.size() - 1), 0) : spline_sections_[i].D;
				if (fitValue < min) // !min.IsEmpty && 
				{
					spline_sections_[i].Y = min + overCorrection * spline_sections_[i].Weight;
					flagModified = true;
				}
				if (fitValue > max) // !max.IsEmpty && 
				{
					spline_sections_[i].Y = max - overCorrection * spline_sections_[i].Weight;
					flagModified = true;
				}
			}

			// Re-initialise spline if modified.
			if (flagModified)
				calculate_Q();
			return flagModified;
		}
	};
}
#endif // cl_tape_examples_impl_spline_cubic_smooth_1d_utils_hpp