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

#ifndef cl_tape_examples_impl_spline_cubic_smooth_1d_hpp
#define cl_tape_examples_impl_spline_cubic_smooth_1d_hpp

#include <cl/tape/tape.hpp>
#include "impl/utils.hpp"
#include "impl/smooth_cubic_spline.hpp"
#include "impl/spline_cubic_smooth_1d_utils.hpp"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

namespace cl
{
	/// <summary>Structure to store spline parameters.</summary>
	class spline_1d_builder
	{
	public:
		// Smoothing parameter.
		tdouble smoothParam_;

		// Fields to bound maximum and minimum allowed spline values.
		tdouble min_;
		tdouble max_;

		// Fields for left and right slopes.
		tdouble slopeLeft_;
		tdouble slopeRight_;

		// Fields for left and right values.
		tdouble valueLeft_;
		tdouble valueRight_;
	public:
		spline_1d_builder() : smoothParam_(0.9), min_(NULL), max_(NULL), slopeLeft_(NULL), slopeRight_(NULL), valueLeft_(NULL), valueRight_(NULL) {}
	}

	///// <summary>Math utility class for ClSplineBase.</summary>
	//class spline_cubic_smooth_1d : smooth_cubic_spline 
	//{
	//private:
	//	// It is true if spline uses in PCA alrogithm
	//	bool _isPca;

	//	///<summary>Apply boundary conditions (add extra sectios).</summary>
	//	void apply_boundary_conditions()
	//	{
	//		// Allow extrapolation.
	//		extrapolate = true;

	//		// Extend spline to match specified slope.
	//		splineMath_.apply_extra_slope(params_.slopeLeft_, params_.slopeRight_);
	//	}
	//protected:

	//	// Underlying spline instance which incorporates all math.
	//	spline_cubic_smooth_1d_utils splineMath_;

	//	/// <summary> Calculates an optimal lambda (rigidity) value by the variance of values_ array.</summary>
	//	tdouble best_rigidity_estimate()
	//	{
	//		if (!is_initialized) return 0;
	//		const double minSmoothParam = 1e-6;
	//		const double maxSmoothParam = 1 - 1e-6;
	//		if (values_.size() < 1) return minSmoothParam;
	//		if (values_.size() < 2) return minSmoothParam;

	//		boost::accumulators::accumulator_set<tdouble, boost::accumulators::features< boost::accumulators::tag::mean, boost::accumulators::tag::variance>, int> sT;
	//		for (auto i = 0; i < values_.size() - 1; ++i)
	//		{
	//			sT(values_[i + 1] - values_[i], boost::accumulators::weight = 1);
	//		}
	//		auto varSt = boost::accumulators::variance(sT);
	//		auto smoothParameter = 1.0 / (varSt * sqrt(values_.size()) + 1);
	//		auto param = NaNConstraint(); // can't Google it :(
	//		if ((smoothParameter > 0.1) && _isPca)
	//			smoothParameter /= sqrt(values_.size());
	//		if ((param.Matches(smoothParameter)) || smoothParameter > maxSmoothParam)
	//			smoothParameter = minSmoothParam;
	//		return smoothParameter;
	//	}

	//	///<summary>Bound spline min and max values.</summary>
	//	void bound(int nIterMax, tdouble overCorrection)
	//	{
	//		//// Set default argument values.
	//		//if (nIterMax.IsEmpty)
	//		//	nIterMax = 100;
	//		//if (overCorrection.IsEmpty)
	//		//	overCorrection = 0.01;

	//		// Iterate while all knots are within bounds or maximum number of iters is reached.
	//		int nIter;
	//		for (nIter = 0; nIter < nIterMax; nIter++)
	//		{
	//			//Console.WriteLine("Iter = " + nIter.ToString() + "  overCorrection = " + overCorrection.ToString());
	//			if (!splineMath_.bound_knots(params_.min_, params_.max_, overCorrection))
	//			{
	//				std::cout << "bound() completed, nIter = " + nIter;
	//				return;
	//			}
	//			splineMath_.solve(params_.smoothParam_);
	//			overCorrection *= 2;
	//		}
	//		std::cout << "bound() failed, reached nIter = " + nIter + (" overCorrection = " + overCorrection.value());
	//	}
	//public:
	//	// Spline parameters.
	//	spline_1d_builder params_;

	//	// Extrapolation flag (provide spline values outside the data).
	//	bool extrapolate;

	//	spline_cubic_smooth_1d(spline_1d_builder parameters) : params_(parameters) {}

	//	/// <summary>Initialise spline with input data.</summary>
	//	void load_data(std::vector<weighted_point> data) override
	//	{
	//		// Check the number of data points
	//		if (data.size() < 4)
	//			throw "Number of data points must be >= 4"; // new ClEx("Number of data points must be >= 4");
	//		// Check dimension of input data
	//		if (data.at(0).array_dimension(data) != 1) // weighted_point -> data.at(0)
	//			throw "The dimension of input data must be 1."; // new ClEx("The dimension of input data must be 1.");
	//		// Deep copy.
	//		data_spline_.resize(data.size());
	//		data_spline_ = data;
	//		// Pre-process input data.
	//		data_spline_ = (weighted_point[])ClAbstractDataProcessor.ProcessDataArray(new ClDataMerge(), data_spline_);
	//		// Set remaining fields.
	//		points_number_ = data_spline_.size();
	//		if (points_number_ < 4)
	//			throw "Minimum 4 points with different x-coordinates are required."; // new ClEx("Minimum 4 points with different x-coordinates are required.");
	//		data.at(0).array_1d_to_XY(data_spline_, points_, values_);
	//		// Initialise underlying spline.
	//		splineMath_.load(data_spline_, params_);
	//		// Update status flags.
	//		is_initialized = true;
	//		is_calculated = false;
	//	}

	//	/// <summary>Returns spline value at a specified point.</summary>
	//	tdouble value_at(tdouble point) override
	//	{
	//		if (!is_calculated)
	//			throw "Spline is not calculated."; // new ClEx("Spline is not calculated.");
	//		const int derivativeOrder = 0;
	//		return splineMath_.value_at(point, derivativeOrder, extrapolate);
	//	}

	//	/// <summary>Returns spline values at input X data.</summary>
	//	std::vector<tdouble> values_at_input() override
	//	{
	//		if (!is_calculated)
	//		throw "Spline is not calculated."; // new ClEx("Spline is not calculated.");

	//		// Array to calculate spline values  at the  reference input X array.
	//		auto indexLast = data_spline_.size() - 1;
	//		std::vector<tdouble> splineAtInput(indexLast + 1);
	//		for (auto i = 0; i < indexLast; i++)
	//			splineAtInput[i] = splineMath_.spline_sections_[i].D;
	//		// Calculate the spline value at the very last right point.
	//		splineAtInput[indexLast] = value_at(splineMath_.spline_sections_[indexLast].X);

	//		return splineAtInput;
	//	}

	//	/// <summary>Calculates spline parameters and extends spline if necessary.</summary>
	//	void calculate() override
	//	{
	//		calculate_without_extension();

	//		// Set boundary conditions if needed.
	//		// TODO. Check whether this can be done BEFORE bound().
	//		if (!params_.slopeLeft_.IsEmpty || !params_.slopeRight_.IsEmpty)
	//			apply_boundary_conditions();

	//		is_calculated = true;
	//	}

	//	/// <summary>Calculates spline parameters without spline extension (applying boundary conditions).</summary>
	//	void calculate_without_extension()
	//	{
	//		if (!is_initialized)
	//			throw "Spline is not initialised."; // new ClEx("Spline is not initialised.");
	//		// Estimate the smoothing parameter
	//		params_.smoothParam_ = best_rigidity_estimate();
	//		// Solve underlying linear system to find spline coefficients.
	//		splineMath_.solve(params_.smoothParam_);

	//		// bound spline if needed.
	//		if (!params_.max_.IsEmpty || !params_.min_.IsEmpty)
	//		{
	//			int nIterMax = 0;
	//			tdouble overCorrection = 0;
	//			bound(nIterMax, overCorrection);
	//		}

	//		is_calculated = true;
	//	}

	//	///// <summary>Return ClCsvMatrix with input data.</summary>
	//	//ClCsvMatrix GetInputData()
	//	//{
	//	//	if (!IsInitialized)
	//	//		throw new ClEx("Spline is not initialised.");

	//	//	// Build matrix with input data.
	//	//	auto matrix = new ClCsvMatrix(points_number_, 2);
	//	//	for (int i = 0; i < points_number_; i++)
	//	//	{
	//	//		matrix[i, 0] = new tdouble(Points[i]).ToVariant();
	//	//		matrix[i, 1] = new tdouble(Values[i]).ToVariant();
	//	//	}
	//	//	return matrix;
	//	//}

	//	///// <summary>Return ClCsvMatrix with fitted data.</summary>
	//	//ClCsvMatrix GetFittedData()
	//	//{
	//	//	if (!IsCalculated)
	//	//		throw new ClEx("Spline is not calculated.");

	//	//	// Build matrix with fitted data.
	//	//	auto matrix = new ClCsvMatrix(splineMath_.DataNumber, 2);
	//	//	for (int i = 0; i < splineMath_.DataNumber; i++)
	//	//	{
	//	//		auto x = Points[i];
	//	//		matrix[i, 0] = new tdouble(x).ToVariant();
	//	//		matrix[i, 1] = new tdouble(ValueAt(x)).ToVariant();
	//	//	}
	//	//	return matrix;
	//	//}

	//	///// <summary>Return ClCsvMatrix with all spline coefficients.</summary>
	//	//ClCsvMatrix GetCoefficients()
	//	//{
	//	//	if (!IsCalculated)
	//	//		throw new ClEx("Spline is not calculated.");
	//	//	return splineMath_.GetCoefficients();
	//	//}
	//};
}
#endif // cl_tape_examples_impl_spline_cubic_smooth_1d_hpp