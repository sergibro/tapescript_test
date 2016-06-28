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
    /// <summary>Math utility class for ClSplineCubicSmoothEndValue1D.</summary>
    public class SplineCubicSmoothEndValue1DUtils : SplineCubicSmooth1DUtils
    {
        protected ClDouble valueLeft_;
        protected ClDouble valueRight_;

        /// <summary> Spline initialisation from X, Y input arrays and boundary conditions (slopes and values at end points)
        /// X must be in ascending order. </summary>
        public override void Load(ClWeightedPoint[] data, Spline1DBuilder parameters)
        {
            // if left or right end value is not set take value from first or last point
            if (parameters.ValueRight.IsEmpty)
                valueRight_ = data[data.Length - 1].Value;
            else valueRight_ = parameters.ValueRight;

            if (parameters.ValueLeft.IsEmpty)
                valueLeft_ = data[0].Value;
            else valueLeft_ = parameters.ValueLeft;

            base.Load(data);

            // Number of unknown b_i.
            equationSize_ = pointsNumber_;
        }

        ///<summary> Initialise right side of equation. </summary>
        protected override void CalculateQ()
        {    
            q_ = new double[pointsNumber_];
            q_[0] = r_[0] * (splineSections_[1].Y - valueLeft_);
            q_[1] = r_[1] * (splineSections_[2].Y - splineSections_[1].Y) - q_[0];
            for (int i = 2; i < pointsNumber_ - 2; ++i)
                q_[i] = 3 * (splineSections_[i + 1].Y - splineSections_[i].Y) / h_[i] - 3 * (splineSections_[i].Y - splineSections_[i - 1].Y) / h_[i - 1];
            q_[pointsNumber_ - 1] = r_[pointsNumber_ - 2] * (splineSections_[pointsNumber_ - 2].Y - valueRight_);
            q_[pointsNumber_ - 2] = r_[pointsNumber_ - 3] * (splineSections_[pointsNumber_ - 3].Y - splineSections_[pointsNumber_ - 2].Y) - q_[pointsNumber_ - 1]; 
        }
        
        ///<summary> Prepare system of linear equations: calculate all non-zero matrix elements.</summary>
        protected override void PrepareLinearSystem(double mu, out double[] u, out double[] v, out double[] w)
        {
            // Diagonal, first and second upper diagonal vectors u, v, w, respectively.
            u = new double[equationSize_];
            v = new double[equationSize_ - 1];
            w = new double[equationSize_ - 2];
            u[0] = mu *  Math.Pow(r_[0], 2) * splineSections_[1].Weight + p_[0];
            v[0] = mu * r_[0] * f_[1] * splineSections_[1].Weight + h_[0];
            w[0] = mu*r_[0]*r_[1]*splineSections_[1].Weight;
            u[1] = mu * (Math.Pow(f_[1], 2) * splineSections_[1].Weight + Math.Pow(r_[1], 2) * splineSections_[2].Weight) + p_[1];
            v[1] = mu * r_[1] * (f_[1] * splineSections_[1].Weight + f_[2] * splineSections_[2].Weight) + h_[1];
            w[1] = mu * r_[1] * r_[2] * splineSections_[2].Weight;
            for (int i = 2; i < equationSize_ - 2; i++)
            {
                u[i] = mu * (Math.Pow(r_[i - 1], 2) * splineSections_[i-1].Weight + Math.Pow(f_[i], 2) * splineSections_[i].Weight
                    + Math.Pow(r_[i], 2) * splineSections_[i + 1].Weight) + p_[i];
                v[i] = mu * r_[i] * (f_[i] * splineSections_[i].Weight + f_[i + 1] * splineSections_[i + 1].Weight) + h_[i];
                w[i] = mu * r_[i] * r_[i + 1] * splineSections_[i + 1].Weight;
            }
            u[equationSize_ - 2] = mu * (Math.Pow(r_[equationSize_ - 3], 2) * splineSections_[equationSize_ - 3].Weight 
                + Math.Pow(f_[equationSize_ - 2], 2) * splineSections_[equationSize_ - 2].Weight) + p_[equationSize_ - 2];
            u[equationSize_ - 1] = mu * Math.Pow(r_[equationSize_ - 2], 2) * splineSections_[equationSize_ - 2].Weight + p_[equationSize_ - 1];
            v[equationSize_ - 2] = mu*r_[equationSize_ - 2]*f_[equationSize_ - 2]*splineSections_[equationSize_ - 2].Weight  + h_[equationSize_ - 2];
        }
        
        
        ///<summary> Calculate all spline coefficients.</summary>
        protected override void FindAllCoefficients(double mu)
        {
            splineSections_[0].D = valueLeft_;
            splineSections_[1].D = splineSections_[1].Y - mu * splineSections_[1].Weight * (r_[0] * q_[0] + f_[1] * q_[1] + r_[1] * q_[2]);
            splineSections_[2].D = splineSections_[2].Y - mu * splineSections_[2].Weight * (r_[1] * q_[1] + f_[2] * q_[2] + r_[2] * q_[3]);
            splineSections_[1].C = (splineSections_[2].D - splineSections_[1].D) / h_[1] - h_[1] * (2 * q_[1] + q_[2]) / 3;
            splineSections_[0].C = 2*h_[0]*(q_[1] + q_[0]) - 2*splineSections_[1].C;            
            splineSections_[1].B = q_[1];
            splineSections_[0].B = q_[0] - 1.5 * splineSections_[0].C / h_[0];
            splineSections_[0].A = (q_[1] - splineSections_[0].B) / (3 * h_[0]); 
            splineSections_[1].A = (q_[2] - q_[1]) / (3 * h_[1]);  
            for (int i = 2; i < pointsNumber_ - 2; ++i)
            {           
                splineSections_[i].B = q_[i];
                splineSections_[i].D = splineSections_[i].Y - mu * splineSections_[i].Weight * (r_[i - 1] * q_[i - 1] + f_[i] * q_[i] + r_[i] * q_[i + 1]);
                splineSections_[i].A = (q_[i + 1] - q_[i]) / (3 * h_[i]);
                splineSections_[i].C = (q_[i] + q_[i - 1]) * h_[i - 1] + splineSections_[i - 1].C;
            }           
            splineSections_[pointsNumber_ - 2].B = q_[pointsNumber_ - 2];
            splineSections_[pointsNumber_ - 2].D = splineSections_[pointsNumber_ - 2].Y - mu * splineSections_[pointsNumber_ - 2].Weight 
                * (r_[pointsNumber_ - 3] * q_[pointsNumber_ - 3] + f_[pointsNumber_ - 2] * q_[pointsNumber_ - 2] + r_[pointsNumber_ - 2] * q_[pointsNumber_ - 1]);
            splineSections_[pointsNumber_ - 1].D = valueRight_;   
            splineSections_[pointsNumber_ - 2].C = (q_[pointsNumber_ - 2] + q_[pointsNumber_ - 3]) * h_[pointsNumber_ - 3] + splineSections_[pointsNumber_ - 3].C;
            splineSections_[pointsNumber_ - 1].C = - 2 * (q_[pointsNumber_ - 1] + q_[pointsNumber_ - 2]) * h_[pointsNumber_ - 2] - 2 * splineSections_[pointsNumber_ - 2].C;
            splineSections_[pointsNumber_ - 1].B = q_[pointsNumber_ - 1] + 1.5 * splineSections_[pointsNumber_ - 1].C / h_[pointsNumber_ - 2];
            splineSections_[pointsNumber_ - 2].A = (splineSections_[pointsNumber_ - 1].B - q_[pointsNumber_ - 2]) / (3 * h_[pointsNumber_ - 2]); ;          
            splineSections_[pointsNumber_ - 1].A = 0;
        }
    }
}
