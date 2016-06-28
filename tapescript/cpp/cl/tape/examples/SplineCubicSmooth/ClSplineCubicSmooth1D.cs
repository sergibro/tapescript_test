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

using System;
using MathNet.Numerics.Statistics;
using NUnit.Framework.Constraints;

namespace Cl
{
    /// <summary>Structure to store spline parameters.</summary>
    public class Spline1DBuilder
    {
        // Smoothing parameter.
        private ClDouble smoothParam_;

        // Fields to bound maximum and minimum allowed spline values.
        private ClDouble min_;
        private ClDouble max_;

        // Fields for left and right slopes.
        private ClDouble slopeLeft_;
        private ClDouble slopeRight_;

        // Fields for left and right values.
        private ClDouble valueLeft_;
        private ClDouble valueRight_;

        public ClDouble SmoothParam { get { return smoothParam_; } set { smoothParam_ = value; } }

        public ClDouble Min { get { return min_; } set { min_ = value; } }

        public ClDouble Max { get { return max_; } set { max_ = value; } }

        public ClDouble SlopeLeft { get { return slopeLeft_; } set { slopeLeft_ = value; } }

        public ClDouble SlopeRight { get { return slopeRight_; } set { slopeRight_ = value; } }

        public ClDouble ValueLeft { get { return valueLeft_; } set { valueLeft_ = value; } }

        public ClDouble ValueRight { get { return valueRight_; } set { valueRight_ = value; } }

        public Spline1DBuilder()
        {
            SmoothParam = 0.9;
            min_ = ClDouble.Empty;
            max_ = ClDouble.Empty;
            slopeLeft_ = ClDouble.Empty;
            slopeRight_ = ClDouble.Empty;
            valueLeft_ = ClDouble.Empty;
            valueRight_ = ClDouble.Empty;
        }
    }

    public class ClSplineCubicSmooth1D : ClAbstractSpline1D
    {
        // Spline parameters.
        public Spline1DBuilder params_;

        // It is true if spline uses in PCA alrogithm
        private bool _isPca;

        public bool IsPCA
        {
            set { _isPca = value; }
        }

        public Spline1DBuilder Params
        {
            get { return params_; }
            set { params_ = value; }
        }

        // Extrapolation flag (provide spline values outside the data).
        public bool Extrapolate { get; set; }

        // Underlying spline instance which incorporates all math.
        protected SplineCubicSmooth1DUtils splineMath_ = new SplineCubicSmooth1DUtils();
        
        /// <summary> Calculates an optimal lambda (rigidity) value by the variance of values_ array.</summary>
        protected double BestRigidityEstimate()
        {
            if (!IsInitialized) return 0;
            const double minSmoothParam = 1e-6;
            const double maxSmoothParam = 1 - 1e-6;
            if (values_.Length < 1) return minSmoothParam;
            var sT = new double[values_.Length];
            if (values_.Length < 2) return minSmoothParam;
            for (var i = 0; i < sT.Length - 1; i++)
                sT[i] = values_[i + 1] - values_[i];
            var varSt = sT.Variance();
            var smoothParameter = 1.0 / (varSt * Math.Pow(sT.Length, 1.0 / 2) + 1);
            var param = new NaNConstraint();
            if ((smoothParameter > 0.1) && _isPca)
                smoothParameter /= Math.Sqrt(sT.Length);
            if ((param.Matches(smoothParameter)) || smoothParameter > maxSmoothParam)
                smoothParameter = minSmoothParam;
            return smoothParameter;
        }

        public ClSplineCubicSmooth1D(Spline1DBuilder parameters)
        {
            params_ = parameters;
        }

        /// <summary>Initialise spline with input data.</summary>
        public override void LoadData(ClWeightedPoint[] data)
        {
            // Check the number of data points
            if (data.Length < 4)
                throw new ClEx("Number of data points must be >= 4");
            // Check dimension of input data
            if (ClWeightedPoint.ArrayDimension(data) != 1)
                throw new ClEx("The dimension of input data must be 1.");
            // Deep copy.
            dataSpline_ = new ClWeightedPoint[data.Length];
            Array.Copy(data, dataSpline_, data.Length);
            // Pre-process input data.
            dataSpline_ = (ClWeightedPoint[])ClAbstractDataProcessor.ProcessDataArray(new ClDataMerge(), dataSpline_);
            // Set remaining fields.
            PointsNumber = dataSpline_.Length;
            if(PointsNumber < 4)
                throw new ClEx("Minimum 4 points with different x-coordinates are required.");
            ClWeightedPoint.Array1DToXY(dataSpline_, out points_, out values_);
            // Initialise underlying spline.
            splineMath_.Load(dataSpline_, params_);
            // Update status flags.
            IsInitialized = true;
            IsCalculated = false;
        }

        /// <summary>Returns spline value at a specified point.</summary>
        public override double ValueAt(double point)
        {
            if (!IsCalculated)
                throw new ClEx("Spline is not calculated.");
            const int derivativeOrder = 0;
            return splineMath_.ValueAt(point, derivativeOrder, Extrapolate);
        }

        /// <summary>Returns spline values at input X data.</summary>
        public override double[] ValuesAtInput 
        {
            get
            {
                if (!IsCalculated)
                    throw new ClEx("Spline is not calculated.");
            
                // Array to calculate spline values  at the  reference input X array.
                var indexLast = dataSpline_.Length - 1;
                var splineAtInput = new double[indexLast + 1];
                for (var i = 0; i < indexLast; i++)
                    splineAtInput[i] = splineMath_.CubicSmoothSections[i].D;
                // Calculate the spline value at the very last right point.
                splineAtInput[indexLast] = ValueAt(splineMath_.CubicSmoothSections[indexLast].X);

                return splineAtInput;
        }
    }

        /// <summary>Calculates spline parameters and extends spline if necessary.</summary>
        public override void Calculate()
        {
            CalculateWithoutExtension();

            // Set boundary conditions if needed.
            // TODO. Check whether this can be done BEFORE Bound().
            if (!params_.SlopeLeft.IsEmpty || !params_.SlopeRight.IsEmpty)
                ApplyBoundaryConditions();

            IsCalculated = true;
        }

        /// <summary>Calculates spline parameters without spline extension (applying boundary conditions).</summary>
        public void CalculateWithoutExtension()
        {
            if (!IsInitialized)
                throw new ClEx("Spline is not initialised.");
            // Estimate the smoothing parameter
            params_.SmoothParam = BestRigidityEstimate();					
            // Solve underlying linear system to find spline coefficients.
            splineMath_.Solve(params_.SmoothParam);

            // Bound spline if needed.
            if (!params_.Max.IsEmpty || !params_.Min.IsEmpty)
            {
                var nIterMax = new ClInt();
                var overCorrection = new ClDouble();
                Bound(nIterMax, overCorrection);
            }

            IsCalculated = true;
        }

        /// <summary>Return ClCsvMatrix with input data.</summary>
        public ClCsvMatrix GetInputData()
        {
            if (!IsInitialized)
                throw new ClEx("Spline is not initialised.");

            // Build matrix with input data.
            var matrix = new ClCsvMatrix(PointsNumber, 2);
            for (ClInt i = 0; i < PointsNumber; i++)
            {
                matrix[i, 0] = new ClDouble(Points[i]).ToVariant();
                matrix[i, 1] = new ClDouble(Values[i]).ToVariant();
            }
            return matrix;
        }

        /// <summary>Return ClCsvMatrix with fitted data.</summary>
        public ClCsvMatrix GetFittedData()
        {
            if (!IsCalculated)
                throw new ClEx("Spline is not calculated.");

            // Build matrix with fitted data.
            var matrix = new ClCsvMatrix(splineMath_.DataNumber, 2);
            for (ClInt i = 0; i < splineMath_.DataNumber; i++)
            {
                var x = Points[i];
                matrix[i, 0] = new ClDouble(x).ToVariant();
                matrix[i, 1] = new ClDouble(ValueAt(x)).ToVariant();
            }
            return matrix;
        }

        /// <summary>Return ClCsvMatrix with all spline coefficients.</summary>
        public ClCsvMatrix GetCoefficients()
        {
            if (!IsCalculated)
                throw new ClEx("Spline is not calculated.");
            return splineMath_.GetCoefficients();
        }

        ///<summary>Apply boundary conditions (add extra sectios).</summary>
        private void ApplyBoundaryConditions()
        {
            // Allow extrapolation.
            Extrapolate = true;

            // Extend spline to match specified slope.
            splineMath_.ApplyExtraSlope(params_.SlopeLeft, params_.SlopeRight);
        }

        ///<summary>Bound spline min and max values.</summary>
        protected void Bound(ClInt nIterMax, ClDouble overCorrection)
        {
            // Set default argument values.
            if (nIterMax.IsEmpty)
                nIterMax = 100;
            if (overCorrection.IsEmpty)
                overCorrection = 0.01;

            // Iterate while all knots are within bounds or maximum number of iters is reached.
            int nIter;
            for (nIter = 0; nIter < nIterMax; nIter++)
            {
                //Console.WriteLine("Iter = " + nIter.ToString() + "  overCorrection = " + overCorrection.ToString());
                if (!splineMath_.BoundKnots(params_.Min, params_.Max, overCorrection))
                {
                    Console.WriteLine("Bound() completed, nIter = " + nIter);
                    return;
                }
                splineMath_.Solve(params_.SmoothParam);
                overCorrection *= 2;
            }
            Console.WriteLine("Bound() failed, reached nIter = " + nIter + " overCorrection = " + overCorrection);
        }
    }
}
