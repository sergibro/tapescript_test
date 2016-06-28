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

namespace Cl
{
    /// <summary>Abstract class for 1D spline.</summary>
    public abstract class ClAbstractSpline1D : IClSpline1D
    {

        protected double[] points_;
        protected double[] values_;
        //Guess function values at input X points.
        protected double[] guessValuesAtInput_;
        // Spline data.
        protected ClWeightedPoint[] dataSpline_;
       
        public bool IsInitialized { get; protected set; }
        public bool IsCalculated { get; protected set; }
        public int PointsNumber { get; protected set; }       
        public double[] Points { get { return points_; } }    
        public double[] Values { get { return values_; } set { throw new ClEx("Use LoadData(ClPoint[]) or LoadData(ClWeightedPoint[])"); } }
        public double[] GuessValuesAtInput { get { return guessValuesAtInput_; } set { guessValuesAtInput_ = value; }}
        public ClPoint[] Data { get { return dataSpline_; } }
        public Spline1DBuilder Params { get; set; }
        public bool Extrapolate { get; set; }

        /// <summary>Initialise spline with unweighted input data.</summary>
        public virtual void LoadData(ClPoint[] data)
        {
            ClWeightedPoint[] dataWeighted = ClWeightedPoint.ToWeightedPointsArray(data);
            LoadData(dataWeighted);
        }

        /// <summary>Initialise spline with weighted input data.</summary>
        public virtual void LoadData(ClWeightedPoint[] data)
        {
            // Check dimension of input data
            if (ClWeightedPoint.ArrayDimension(data) != 1)
                throw new ClEx("The dimension of input data must be 1.");
            // Deep copy.
            dataSpline_ = new ClWeightedPoint[data.Length];
            Array.Copy(data, dataSpline_, data.Length);
            // Set remaining fields.
            PointsNumber = dataSpline_.Length;
            ClWeightedPoint.Array1DToXY(dataSpline_, out points_, out values_);
        }

        /// <summary>Calculates spline parameters.</summary>
        public abstract void Calculate();

        /// <summary>Returns spline value at a specified point.</summary>
        public abstract double ValueAt(double point);

        /// <summary>Returns spline values at  all input X points.</summary>
        public virtual  double[] ValuesAtInput {
            get
            {          
                double[] splineAtInput = new double[PointsNumber];
                for (int i = 0; i < PointsNumber; i++)
                    splineAtInput[i] = ValueAt(dataSpline_[i].Value);
                return splineAtInput;
            } 
        }


        /// <summary>Returns spline values on grid with a specified number of nodes.</summary>
        public ClCsvMatrix ValueOnGrid(int nodeCount)
        {
            ClDouble knotStep = (Points[PointsNumber - 1] - Points[0]) / (nodeCount - 1);
            return ValueOnGrid(knotStep);
        }

        /// <summary>Returns spline values on the grid with a specified step.</summary>
        public ClCsvMatrix ValueOnGrid(double nodeStep)
        {
            ClInt n = Convert.ToInt32(Math.Ceiling((Points[PointsNumber - 1] - Points[0]) / nodeStep));
            ClCsvMatrix matrix = new ClCsvMatrix(n, 2);
            for (ClInt i = 0; i < n; i++)
            {
                ClDouble x = Points[0] + i * nodeStep;
                matrix[i, 0] = x.ToVariant();
                matrix[i, 1] = ValueAt(x);
            }
            return matrix;
        }
    }
}