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
    /// <summary>Math utility class for ClSplineCubicSmoothSlopeValue1D.</summary>
    public class SplineCubicSmoothSlopeValue1DUtils : SplineCubicSmoothEndValue1DUtils
    {
        private ClDouble slopeLeft_;
        private ClDouble slopeRight_;

        /// <summary> Spline initialisation from X, Y input arrays and boundary conditions (values at end points)
        /// X must be in ascending order. </summary>
        public override void Load(ClWeightedPoint[] data, Spline1DBuilder parameters)
        {
            // if left or right slope is not set calculate it from first or last two points
            if (parameters.SlopeRight.IsEmpty)
                slopeRight_ = (data[data.Length - 1].Value - data[data.Length - 2].Value) /
                              (data[data.Length - 1].X[0] - data[data.Length - 2].X[0]);
            else slopeRight_ = parameters.SlopeRight;

            if (parameters.SlopeLeft.IsEmpty)
                slopeLeft_ = (data[1].Value - data[0].Value) / (data[1].X[0] - data[0].X[0]);
            else slopeLeft_ = parameters.SlopeLeft;

            base.Load(data, parameters);
        }

        ///<summary> Initialise right side of equation. </summary>
        protected override void CalculateQ()
        {    
            base.CalculateQ();
            q_[0] += -3*slopeLeft_;
            q_[pointsNumber_ - 1] += 3*slopeRight_;
        }
              
        ///<summary> Calculate all spline coefficients.</summary>
        protected override void FindAllCoefficients(double mu)
        {
            splineSections_[0].B = q_[0];
            splineSections_[0].D = valueLeft_;
            splineSections_[0].A = (q_[1] - q_[0]) / (3 * h_[0]);
            splineSections_[0].C = slopeLeft_;
            for (int i = 1; i < pointsNumber_ - 1; ++i)
            {
                splineSections_[i].B = q_[i];
                splineSections_[i].D = splineSections_[i].Y - mu * splineSections_[i].Weight * (r_[i - 1] * q_[i - 1] + f_[i] * q_[i] + r_[i] * q_[i + 1]);
                splineSections_[i].A = (q_[i + 1] - q_[i]) / (3 * h_[i]);
                splineSections_[i].C = (q_[i] + q_[i - 1]) * h_[i - 1] + splineSections_[i - 1].C;
            }
            splineSections_[pointsNumber_ - 1].B = q_[pointsNumber_ - 1];
            splineSections_[pointsNumber_ - 1].D = valueRight_;
            splineSections_[pointsNumber_ - 1].A = 0;
            splineSections_[pointsNumber_ - 1].C = slopeRight_;
        }
    }
}
