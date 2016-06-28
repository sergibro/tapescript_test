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

#ifndef cl_tape_examples_impl_abstract_spline_1d_hpp
#define cl_tape_examples_impl_abstract_spline_1d_hpp

#include <cl/tape/tape.hpp>
#include "impl/utils.hpp"

namespace cl
{
	/// <summary>Abstract class for 1D spline.</summary>
	class abstract_spline_1d
	{
	protected:
		std::vector<tdouble> points_;
		std::vector<tdouble> values_;
		//Guess function values at input X points.
		std::vector<tdouble> guess_values_at_input_;
		// Spline data.
		std::vector<weighted_point> data_spline_;

		bool is_initialized;
		bool is_calculated;
		int points_number;
	public:
		spline_1d_builder params;
		bool extrapolate;
		bool get_is_initialized(){ return is_initialized; }
		bool get_is_calculated(){ return is_calculated; }
		int get_points_number(){ return points_number; }
		std::vector<tdouble> Points(){ return points_; }
		std::vector<tdouble> get_values(){ return values_; }
		//void set_values(std::vector<tdouble>) { throw new ClEx("Use load_data(ClPoint[]) or load_data(ClWeightedPoint[])"); }
		std::vector<tdouble> get_guess_values_at_input(){ return guess_values_at_input_; }
		void set_guess_values_at_input(std::vector<tdouble> val) { guess_values_at_input_ = val; }
		std::vector<point> get_data() { return data_spline_; }

		/// <summary>Initialise spline with unweighted input data.</summary>
		virtual void load_data(std::vector<point> data)
		{
			std::vector<weighted_point> data_weighted = data;
			load_data(data_weighted);
		}

		/// <summary>Initialise spline with weighted input data.</summary>
		virtual void load_data(std::vector<weighted_point> data)
		{
			// Check dimension of input data
			if (data.at(0).array_dimension(data) != 1) // weighted_point -> data.at(0)
				throw "The dimension of input data must be 1."; // ClEx("The dimension of input data must be 1.");
			// Deep copy.
			data_spline_.resize(data.size());
			data_spline_ = data;
			// Set remaining fields.
			points_number = data_spline_.size();
			data.at(0).array_1d_to_XY(data_spline_, points_, values_); // weighted_point -> data.at(0)
		}

		/// <summary>Calculates spline parameters.</summary>
		virtual void calculate()=0;

		/// <summary>Returns spline value at a specified point.</summary>
		virtual tdouble value_at(tdouble point)=0;

		/// <summary>Returns spline values at  all input X points.</summary>
		virtual  std::vector<tdouble> values_at_input(){
			std::vector<tdouble> spline_at_input(points_number);
			for (int i = 0; i < points_number; i++)
				spline_at_input[i] = value_at(data_spline_[i].y);
			return spline_at_input;
		}

		///// <summary>Returns spline values on grid with a specified number of nodes.</summary>
		//ClCsvMatrix value_on_grid(int node_count)
		//{
		//	ClDouble knot_step = (Points()[points_number - 1] - Points()[0]) / (node_count - 1);
		//	return value_on_grid(knot_step);
		//}

		///// <summary>Returns spline values on the grid with a specified step.</summary>
		//ClCsvMatrix value_on_grid(double node_step)
		//{
		//	ClInt n = Convert.ToInt32(Math.Ceiling((Points()[points_number - 1] - Points()[0]) / node_step));
		//	ClCsvMatrix matrix = new ClCsvMatrix(n, 2);
		//	for (ClInt i = 0; i < n; i++)
		//	{
		//		ClDouble x = Points()[0] + i * node_step;
		//		matrix[i, 0] = x.ToVariant();
		//		matrix[i, 1] = value_at(x);
		//	}
		//	return matrix;
		//}

	};
}
#endif // cl_tape_examples_impl_abstract_spline_1d_hpp