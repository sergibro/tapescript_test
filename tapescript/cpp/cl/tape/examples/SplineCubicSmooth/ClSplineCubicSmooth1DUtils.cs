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
using System.Collections.Generic;
using System.Linq;

namespace Cl
{
    /// <summary>Math utility class for ClSplineBase.</summary>
    public class SplineCubicSmooth1DUtils
    {
        // Number of data point and spline sections.
        protected int pointsNumber_;
        public int DataNumber { get { return pointsNumber_; } set { pointsNumber_ = value; } }

        // Array of spline sections.
        protected ClSplineCubicSection[] splineSections_;
        public ClSplineCubicSection[] CubicSmoothSections { get { return splineSections_; } set { splineSections_ = value; } }

        // Helper variables needed for fast and compact calculations.
        // (their one-letter names exactly follow Pollock's article).
        protected double[] p_;
        protected double[] h_;
        protected double[] r_;
        protected double[] f_;
        protected double[] q_;

        // Dimension of the system of linear equations (number of unknown b_i).
        protected int equationSize_;

        /// <summary> Spline initialisation from X, Y input arrays.
        /// u is array of uncertainties.
        /// X must be in ascending order. </summary>
        virtual public void Load(ClWeightedPoint[] data, Spline1DBuilder  parameters = null)
        {
            // Number of points.
            pointsNumber_ = data.Length;

            // Number of unknown b_i.
            equationSize_ = pointsNumber_ - 2;

            // Initialise spline sections (input points are (x_0, y_0), ..., (x_n, y_n)).
            splineSections_ = new ClSplineCubicSection[pointsNumber_];
            for (int i = 0; i < pointsNumber_; ++i)
            {
                splineSections_[i] = new ClSplineCubicSection {X = data[i].X[0], Y = data[i].Value, Weight = 1};
                // TODO Fix NaN problem with weigths.
                //splineSections_[i].Weight = data[i].Weight;
            }

            // Calculate spline helper variables.
            CalculateHelperVariables();
        }

        /// <summary>Returns knot section coordinate.</summary>
        public double GetKnotValue(int section)
        {
            return splineSections_[section].X;
        }

        /// <summary>Returns data value in section.</summary>
        public double GetDataValue(int section)
        {
            return splineSections_[section].Y;
        }

        /// <summary>Get spline derivative at point x with possible extrapolation.</summary>
        public double ValueAt(double x, int derivative, bool extrapolate)
        {
            // If spline can be extrapolated, check if requested x is outside its sections.
            if (extrapolate)
            {
                if (x < splineSections_[0].X)
                    return GetSplineValueInSection(0, x - splineSections_[0].X, derivative);
                else if (x > splineSections_[splineSections_.Length - 1].X)
                    return GetSplineValueInSection(splineSections_.Length - 2, x - splineSections_[splineSections_.Length - 2].X, derivative);
            }

            // Check for requested x within spline.
            // TODO Use binary search here.
            for (int j = 0; j < splineSections_.Length - 1; ++j)
            {
                if (x >= splineSections_[j].X && x <= splineSections_[j + 1].X)
                    return GetSplineValueInSection(j, x - splineSections_[j].X, derivative);
            }

            // If spline section for requested x was not found, return empty.
            throw new ClEx("Value not found.");
        }

        ///<summary> Return spline value in the specified section (j) at the specified internal point (x).
        /// Optionaly specify derivative (0 by default, i.e. spline value).</summary>
        public double GetSplineValueInSection(int j, double x, int derivative)
        {
            return splineSections_[j].ValueAt(x, derivative);
        }

        /// <summary>  Return ClCsvMatrix with all spline coefficients.  </summary>
        public ClCsvMatrix GetCoefficients()
        {
            ClCsvMatrix matrix = new ClCsvMatrix(splineSections_.Length - 1, 6);
            for (int i = 0; i < splineSections_.Length - 1; i++)
            {
                matrix[i, 0] = new ClDouble(splineSections_[i].X).ToVariant();
                matrix[i, 1] = new ClDouble(splineSections_[i + 1].X).ToVariant();
                matrix[i, 2] = new ClDouble(splineSections_[i].A).ToVariant();
                matrix[i, 3] = new ClDouble(splineSections_[i].B).ToVariant();
                matrix[i, 4] = new ClDouble(splineSections_[i].C).ToVariant();
                matrix[i, 5] = new ClDouble(splineSections_[i].D).ToVariant();
            }
            return matrix;
        }

        ///<summary> Calculate spline coefficients.</summary>
        public void Solve(double lambda)
        {
            // Bound smoothing parameter value to avoid problems with NaN.
            const double lambdaMin = 1e-6;
            const double lambdaMax = 1 - 1e-6;
            if (lambda < lambdaMin)
                lambda = lambdaMin;
            else if (lambda > lambdaMax)
                lambda = lambdaMax;

            // Internal smoothing parameter: transform lambda -> mu = 2 * (1 - lambda) / (3 * lambda).
            double mu = 2 * (1 - lambda) / (3 * lambda);

            // Prepare system of linear equations: calculate diagonal, first and second upper diagonal vectors u, v, w, respectively.
            double[] u, v, w;
            PrepareLinearSystem(mu, out u, out v, out w);

            // Find b-coefficients in splines in q_ array.
            FindBCoefficients(ref u, ref v, ref w);

            // Find the rest of spline parameters.
            FindAllCoefficients(mu);
        }

        ///<summary>Apply boundary conditions (add extra sectios) 
        /// to match specified slopes slope1 and slope2.</summary>
        public void ApplyExtraSlope(double slope1, double slope2)
        {
            // Create list of existing spline sections.
            List<ClSplineCubicSection> listSection = new List<ClSplineCubicSection>(splineSections_);

            // Add extra sections at the beginning.
            double extraStep = h_.First();
            ClSplineCubicSection[] extraCubicSmoothSectionBegin = listSection[0].Extend(slope1, extraStep);
            for (int i = 2; i >= 0; i--)
                listSection.Insert(0, extraCubicSmoothSectionBegin[i]);

            // Add extra sections at the end (this is less trivial).
            // First calculate negative extra step.
            extraStep = -h_.Last();
            ClSplineCubicSection cubicSmoothSectionLast = listSection[listSection.Count - 2];
            // Re-reference last spline section.
            cubicSmoothSectionLast.ChangeStartPoint(listSection.Last().X, listSection.Last().Y);
            // Obtain needed extra sections.
            ClSplineCubicSection[] extraCubicSmoothSectionEnd = cubicSmoothSectionLast.Extend(slope2, extraStep);
            // And re-reference them all back.
            ClSplineCubicSection cubicSmoothSectionTail = new ClSplineCubicSection
            {
                X = extraCubicSmoothSectionEnd[0].X,
                Y = extraCubicSmoothSectionEnd[0].Y
            };
            extraCubicSmoothSectionEnd[0].ChangeStartPoint(extraCubicSmoothSectionEnd[1].X, extraCubicSmoothSectionEnd[1].Y);
            extraCubicSmoothSectionEnd[1].ChangeStartPoint(extraCubicSmoothSectionEnd[2].X, extraCubicSmoothSectionEnd[2].Y);
            extraCubicSmoothSectionEnd[2].ChangeStartPoint(listSection[listSection.Count - 1].X, listSection[listSection.Count - 1].Y);
            // Finally, remove arteficial last section and add new ones.
            listSection.RemoveAt(listSection.Count - 1);
            for (int i = 2; i >= 0; i--)
                listSection.Add(extraCubicSmoothSectionEnd[i]);
            listSection.Add(cubicSmoothSectionTail);

            // Update field splineSections_.
            splineSections_ = listSection.ToArray();
        }

        ///<summary> Bound spline knots.</summary>
        public bool BoundKnots(ClDouble min, ClDouble max, ClDouble overCorrection)
        {
            // Set default argument value.
            if (overCorrection.IsEmpty)
                overCorrection = 0;

            // Check knots, correct data if needed.
            bool flagModified = false;
            for (int i = 0; i < pointsNumber_; i++)
            {
                double fitValue = (i == pointsNumber_ - 1) ? GetSplineValueInSection(i - 1, h_.Last(), 0) : splineSections_[i].D;
                if (!min.IsEmpty && fitValue < min)
                {
                    splineSections_[i].Y = min + overCorrection * splineSections_[i].Weight;
                    flagModified = true;
                }
                if (!max.IsEmpty && fitValue > max)
                {
                    splineSections_[i].Y = max - overCorrection * splineSections_[i].Weight;
                    flagModified = true;
                }
            }

            // Re-initialise spline if modified.
            if (flagModified)
                CalculateQ();
            return flagModified;
        }

        ///<summary> Initialise spline helper variables. </summary>
        protected void CalculateHelperVariables()
        {
            p_ = new double[pointsNumber_];
            h_ = new double[pointsNumber_ - 1];
            r_ = new double[pointsNumber_ - 1];
            f_ = new double[pointsNumber_];
            h_[0] = splineSections_[1].X - splineSections_[0].X;
            r_[0] = 3 / h_[0];
            p_[0] = 2 * h_[0];
            f_[0] = -r_[0];
            for (int i = 1; i < pointsNumber_ - 1; ++i)
            {
                h_[i] = splineSections_[i + 1].X - splineSections_[i].X;
                r_[i] = 3 / h_[i];
                p_[i] = 2 * (h_[i - 1] + h_[i]);
                f_[i] = -(r_[i - 1] + r_[i]);
            }
            p_[pointsNumber_ - 1] = 2 * h_[pointsNumber_ - 2];
            f_[pointsNumber_ - 1] = -r_[pointsNumber_ - 2];
            CalculateQ();
        }

        ///<summary> Initialise right side of equation. </summary>
        virtual protected void CalculateQ()
        {
            q_ = new double[pointsNumber_ - 2];
            for (int i = 1; i < pointsNumber_ - 1; ++i)
                q_[i - 1] = 3 * (splineSections_[i + 1].Y - splineSections_[i].Y) / h_[i] - 3 * (splineSections_[i].Y - splineSections_[i - 1].Y) / h_[i - 1];
        }

        ///<summary> Prepare system of linear equations: calculate all non-zero matrix elements.</summary>
        virtual protected void PrepareLinearSystem(double mu, out double[] u, out double[] v, out double[] w)
        {
            // Diagonal, first and second upper diagonal vectors u, v, w, respectively.
            u = new double[equationSize_];
            v = new double[equationSize_ - 1];
            w = new double[equationSize_ - 2];
            for (int i = 0; i < equationSize_ - 2; i++)
            {
                u[i] = mu * (Math.Pow(r_[i], 2) * splineSections_[i].Weight + Math.Pow(f_[i + 1], 2) * splineSections_[i + 1].Weight
                    + Math.Pow(r_[i + 1], 2) * splineSections_[i + 2].Weight) + p_[i + 1];
                v[i] = mu * (r_[i + 1] * (f_[i + 1] * splineSections_[i + 1].Weight + f_[i + 2] * splineSections_[i + 2].Weight)) + h_[i + 1];
                w[i] = mu * r_[i + 1] * r_[i + 2] * splineSections_[i + 2].Weight;
            }
            u[equationSize_ - 2] = mu * (Math.Pow(r_[equationSize_ - 2], 2) * splineSections_[equationSize_ - 2].Weight +
                Math.Pow(f_[equationSize_ - 1], 2) * splineSections_[equationSize_ - 1].Weight +
                Math.Pow(r_[equationSize_ - 1], 2) * splineSections_[equationSize_].Weight) + p_[equationSize_ - 1];
            u[equationSize_ - 1] = mu * (Math.Pow(r_[equationSize_ - 1], 2) * splineSections_[equationSize_ - 1].Weight +
                Math.Pow(f_[equationSize_], 2) * splineSections_[equationSize_].Weight +
                Math.Pow(r_[equationSize_], 2) * splineSections_[equationSize_ + 1].Weight) + p_[equationSize_];
            v[equationSize_ - 2] = mu * r_[equationSize_ - 1] * (f_[equationSize_ - 1] * splineSections_[equationSize_ - 1].Weight +
                f_[equationSize_] * splineSections_[equationSize_].Weight) + h_[equationSize_ - 1];
        }

        ///<summary> Calculate spline B coefficients 
        /// (solution of linear system of equations).</summary>
        protected void FindBCoefficients(ref double[] u, ref double[] v, ref double[] w)
        {
            // Factorization procedure.
            w[0]/=u[0];
            v[0]/=u[0];
            u[1]-= u[0] * Math.Pow(v[0], 2);
            v[1]=(v[1] - u[0] * v[0] * w[0]) / u[1];
            w[1]/=u[1];
            for (int i = 3; i <= equationSize_ - 2; i++)
            {
                u[i - 1]-=(u[i - 3] * Math.Pow(w[i - 3], 2)+u[i - 2] * Math.Pow(v[i - 2], 2));
                v[i - 1] = (v[i - 1] - u[i - 2] * v[i - 2] * w[i - 2]) / u[i - 1];
                w[i - 1]/=u[i - 1];
            }
            u[equationSize_ - 2] -= (u[equationSize_ - 4] * Math.Pow(w[equationSize_ - 4], 2) +
                u[equationSize_ - 3] * Math.Pow(v[equationSize_ - 3], 2));
            v[equationSize_ - 2] = (v[equationSize_ - 2] - u[equationSize_ - 3] * v[equationSize_ - 3] * w[equationSize_ - 3]) / u[equationSize_ - 2];
            u[equationSize_ - 1] -= (u[equationSize_ - 3] * Math.Pow(w[equationSize_ - 3], 2) +
                u[equationSize_ - 2] * Math.Pow(v[equationSize_ - 2], 2));

            // Forward substitution.
            q_[1]-=v[0] * q_[0];
            for (int i = 2; i < equationSize_; ++i)
                q_[i]-=(v[i - 1] * q_[i - 1] + w[i - 2] * q_[i - 2]);
            for (int i = 0; i < equationSize_; ++i)
                q_[i]/= u[i];

            // Backward substitution.
            q_[equationSize_ - 2] = q_[equationSize_ - 2] - v[equationSize_ - 2] * q_[equationSize_ - 1];
            for (int i = equationSize_ - 3; i >= 0; --i)
                q_[i]-= (v[i] * q_[i + 1]+ w[i] * q_[i + 2]);
        }

        ///<summary> Calculate all spline coefficients.</summary>
        virtual protected void FindAllCoefficients(double mu)
        {
            splineSections_[0].D = splineSections_[0].Y - mu * r_[0] * q_[0] * splineSections_[0].Weight;
            splineSections_[1].D = splineSections_[1].Y - mu * splineSections_[0].Weight * (f_[1] * q_[0] + r_[1] * q_[1]);
            splineSections_[0].A = q_[0] / (3 * h_[0]);
            splineSections_[0].B = 0.0;
            splineSections_[0].C = (splineSections_[1].D - splineSections_[0].D) / h_[0] - q_[0] * h_[0] / 3;
            splineSections_[1].A = (q_[1] - q_[0]) / (3 * h_[1]);
            splineSections_[1].B = q_[0];
            splineSections_[1].C = q_[0] * h_[0] + splineSections_[0].C;
            for (int i = 2; i < pointsNumber_ - 2; ++i)
            {
                splineSections_[i].A = (q_[i] - q_[i - 1]) / (3 * h_[i]);
                splineSections_[i].B = q_[i - 1];
                splineSections_[i].C = (q_[i - 1] + q_[i - 2]) * h_[i - 1] + splineSections_[i - 1].C;
                splineSections_[i].D = splineSections_[i].Y - mu * splineSections_[i].Weight * (r_[i - 1] * q_[i - 2] + f_[i] * q_[i - 1] + r_[i] * q_[i]);
            }
            splineSections_[pointsNumber_ - 2].A = (0 - q_[pointsNumber_ - 3]) / (3 * h_[pointsNumber_ - 2]);
            splineSections_[pointsNumber_ - 2].B = q_[pointsNumber_ - 3];
            splineSections_[pointsNumber_ - 2].C = (q_[pointsNumber_ - 3] + q_[pointsNumber_ - 4]) * h_[pointsNumber_ - 3] + splineSections_[pointsNumber_ - 3].C;
            splineSections_[pointsNumber_ - 2].D = splineSections_[pointsNumber_ - 2].Y - mu * splineSections_[pointsNumber_ - 2].Weight 
                * (r_[pointsNumber_ - 3] * q_[pointsNumber_ - 4] + f_[pointsNumber_ - 2] * q_[pointsNumber_ - 3]);
        }
    }
}
