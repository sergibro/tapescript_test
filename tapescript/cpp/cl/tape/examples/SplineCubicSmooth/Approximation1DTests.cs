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
using System.Linq;
using NUnit.Framework;

namespace Cl.Engine.Finance.Approximation.Tests
{
    public class Approximation1DTests : ClUnitTests
    {
        // Additional steering parameters.
        private const int StartSeed = 0;
        private const bool PrintPerformance = true;
        // Control production of output.
        const bool PlotResults = true;
        const bool ConsoleOutput = true;

        /// <summary>Test polynomial regression (class ClPolynomialLsqAmcBuilder).
        /// Input *.csv files: {dataX, dataY, refY},
        ///       where dataX and dataY are input data values, refY is reference.
        /// Output files:
        ///    - *SplineData.csv {dataX, dataY} (input data)
        ///    - *SplineValues.csv {dataX, dataY, refY, fitY},
        ///      where fitY is spline value (i.e. approximated value);
        /// Note that {dataX, dataY} in input and output files may differ because of 
        /// sorting and merging data points with the same x.</summary>
        
        [Test]
        public void AMCRegressionTest()
        {
            ClString name = "AMCRegressionTest";
            BeginTest(name);

            // Testing data: Call, Put, Risk Reversal, Straddle, Strangle, Bull, Bear, Butterfly and Barrier (Down-and-in call, Down-and-out call, Up-and-in call, Up-and-out call ) Options.

            ApproximationTestsHelper test = new ApproximationTestsHelper(Log, Context, name, ConsoleOutput);
            foreach (ClString testName in test.GetInputFileNames(GetType().FullName))
            {
                test.StartTest(testName);
                
                // Read input data.
                ClString fileName = testName + ".csv";
                ClPoint[] data, reference;
                test.ReadCSVInputData(fileName, out data, out reference);
                // initial data for polynomial regression
                double[] dataX = new double[data.Length];
                for (int pos = 0; pos < data.Length; pos++)
                    dataX[pos] = data[pos].X[0];
                ClDoubleArray explanatoryVariables = new ClDoubleArray(dataX);
                double[] dataV = new double[data.Length];
                for (int pos = 0; pos < data.Length; pos++)
                    dataV[pos] = data[pos].Value;
                ClDoubleArray exposureDiscounted = new ClDoubleArray(dataV);
                // Initialise polynomial regression and calculate it
                var amcBuilder = new ClPolynomialLsqAmcBuilder
                {
                    Sdf = new ClDoubleArray(data.Length, 1.0),
                    ExplanatoryVariables = explanatoryVariables,
                    PolynomialPower = 3
                };
                IClAmc amc = amcBuilder.Build<IClAmc>();

                test.CalculationTime.Restart();
                ClDoubleArray result = new ClDoubleArray(amc.Evaluate(exposureDiscounted));
                test.CalculationTime.Stop();

                ClPoint[] valuesAt = new ClPoint[data.Length];
                for (int pos = 0; pos < data.Length; pos++)
                    valuesAt[pos] = new ClPoint(data[pos].X, result[pos]);

                // Initialise spline and calculate spline coefficients.
                var spline = new ClSplineCubicSmooth1D(new Spline1DBuilder());
                spline.LoadData(data);
                spline.Calculate();
                
                // Store results and complete test.
                if (PlotResults)
                    test.StoreSpline1D(spline, valuesAt, data, reference);
                test.Print("Completed " + name + " for " + testName);
            }
            test.FinishTest(PrintPerformance);
        }

        /// <summary>Test smoothed cubic B-splines (class Cl.ClSplineCubicSmooth1D).
        /// Input *.csv files: {dataX, dataY, refY},
        ///       where dataX and dataY are input data values, refY is reference.
        /// Output files:
        ///    - *SplineData.csv {dataX, dataY} (input data)
        ///    - *SplineValues.csv {dataX, dataY, refY, fitY},
        ///      where fitY is spline value (i.e. approximated value);
        ///    - *SplineCoefficients.csv {dataXleft, dataXright, a, b, c, d},
        ///      where a, b, c, d are spline coefficients.
        /// Note that {dataX, dataY} in input and output files may differ because of 
        /// sorting and merging data points with the same x.</summary>
        
       [Test]
        public void SplineTest()
        {
            ClString name = "SplineTest";
            BeginTest(name);

            // Testing data: Call, Put, Risk Reversal, Straddle, Strangle, Bull, Bear, Butterfly and Barrier (Down-and-in call, Down-and-out call, Up-and-in call, Up-and-out call ) Options.
            
            ApproximationTestsHelper test = new ApproximationTestsHelper(Log, Context, name, ConsoleOutput);
            foreach (ClString testName in test.GetInputFileNames(GetType().FullName))
            {
                test.StartTest(testName);

                // Read input data.
                ClString fileName = testName + ".csv";
                ClPoint[] data, reference;
                test.ReadCSVInputData(fileName, out data, out reference);

                ClSplineCubicSmooth1D spline = new ClSplineCubicSmooth1D(new Spline1DBuilder());
                spline.LoadData(data);

                test.CalculationTime.Restart();
                spline.Calculate();
                test.CalculationTime.Stop();

                // Store results and complete test.
                test.Print(String.Format("Optimal lambda: {0}", spline.Params.SmoothParam));
                if (PlotResults)
                    test.StoreSpline1D(spline, reference);
                test.Print("Completed " + name + " for " + testName);
            }
            test.FinishTest(PrintPerformance);
        }

        /// <summary>Test bounded smoothed cubic B-splines (class Cl.ClSplineCubicSmooth1D).
        /// Input:
        ///    - *.csv files: {dataX, dataY, refY},
        ///       where dataX and dataY are input data values, refY is reference.
        ///    - *Bounds.csv {min max},
        ///       where min and max are required minimum and maximum boundaries of the spline, respectively.
        /// Output files:
        ///    - *SplineData.csv {dataX, dataY} (input data)
        ///    - *SplineValues.csv {dataX, dataY, refY, fitY},
        ///      where fitY is spline value (i.e. approximated value);
        ///    - *SplineCoefficients.csv {dataXleft, dataXright, a, b, c, d},
        ///      where a, b, c, d are spline coefficients.
        /// Note that {dataX, dataY} in input and output files may differ because of 
        /// sorting and merging data points with the same x.</summary>
        [Test]
        public void SplineBoundTest()
        {
            ClString name = "SplineBoundTest";
            BeginTest(name);

            // Testing data: CallGeomBrownian, PayoffButterfly, PayoffCallOption.

            var test = new ApproximationTestsHelper(Log, Context, name, ConsoleOutput);
            foreach (ClString testName in test.GetInputFileNames(GetType().FullName, "Bounds"))
            {
                test.StartTest(testName);

                // Read input data.
                ClString fileName = testName + ".csv";
                ClPoint[] data, reference;
                test.ReadCSVInputData(fileName, out data, out reference);

                // Read bounds.
                ClDouble boundMin;
                ClDouble boundMax;
                const bool fromFile = false;
                // TODO Switch to true to load bounds from file
                if (fromFile)
                {
                    fileName = testName + "Bounds.csv";
                    ClCsvMatrix matrix = new ClCsvMatrix(Context.AsFileContext.LoadText(fileName));
                    boundMin = matrix[0, 0].ToDouble();
                    boundMax = matrix[0, 1].ToDouble();
                }
                else
                {
                    var refer = new double[reference.Length];
                    for (var i = 0; i < refer.Length; i++)
                        refer[i] = reference[i].Value;
                    boundMin = refer.Min();
                    boundMax = refer.Max();
                }
                
                // Initialise spline and calculate spline coefficients.
                var splineParams = new Spline1DBuilder {Min = boundMin, Max = boundMax};

                var spline = new ClSplineCubicSmooth1D(splineParams);
                spline.LoadData(data);
                test.CalculationTime.Restart();
                spline.Calculate();
                test.CalculationTime.Stop();

                // Store results and complete test.
                if (PlotResults)
                    test.StoreSpline1D(spline, reference);
                test.Print("Completed " + name + " for " + testName);
            }
            test.FinishTest(PrintPerformance);
        }

        /// <summary>Test smoothed cubic B-splines with extra slope (class Cl.ClSplineCubicSmooth1D).
        /// Input:
        ///    - *.csv files: {dataX, dataY, refY},
        ///       where dataX and dataY are input data values, refY is reference.
        ///    - *Slope.csv {left right},
        ///       where left and right are required left and right slopes of the spline, respectively.
        /// Output files:
        ///    - *SplineData.csv {dataX, dataY} (input data)
        ///    - *SplineValues.csv {dataX, dataY, refY, fitY},
        ///      where fitY is spline value (i.e. approximated value);
        ///    - *SplineCoefficients.csv {dataXleft, dataXright, a, b, c, d},
        ///      where a, b, c, d are spline coefficients.
        /// Note that {dataX, dataY} in input and output files may differ because of 
        /// sorting and merging data points with the same x.</summary>
        [Test]
        public void SplineSlopeTest()
        {
            ClString name = "SplineSlopeTest";
            BeginTest(name);

            // Testing data: CallGeomBrownian; Linear, Parabolic, Exponential Data; PayoffButterfly, PayoffCallOption.

            var test = new ApproximationTestsHelper(Log, Context, name, ConsoleOutput);
            foreach (ClString testName in test.GetInputFileNames(GetType().FullName, "Slope"))
            {
                test.StartTest(testName);

                // Read input data.
                ClString fileName = testName + ".csv";
                ClPoint[] data, reference;
                test.ReadCSVInputData(fileName, out data, out reference);

                // Read slope.
                fileName = testName + "Slope.csv";
                var matrix = new ClCsvMatrix(Context.AsFileContext.LoadText(fileName));
                var slopeLeft = matrix[0, 0].ToDouble();
                var slopeRight = matrix[0, 1].ToDouble();

                // Initialise spline and calculate spline coefficients.
                var spline = new ClSplineCubicSmooth1D(new Spline1DBuilder())
                {
                    Params =
                    {
                        SlopeLeft = slopeLeft,
                        SlopeRight = slopeRight
                    }
                };
                spline.LoadData(data);
                test.CalculationTime.Restart();
                spline.Calculate();
                test.CalculationTime.Stop();
                // Store results and complete test.
                test.Print(String.Format("Optimal lambda: {0}", spline.Params.SmoothParam));
                if (PlotResults)
                    test.StoreSpline1D(spline, reference);
                test.Print("Completed " + name + " for " + testName);
            }
            test.FinishTest(PrintPerformance);
        }

        /// <summary>Test smoothed cubic B-splines with optimal smoothing parameter 
        /// (class Cl.ClSplineCubicOptimalSmooth1D).
        /// Input *.csv files: {dataX, dataY, refY},
        ///       where dataX and dataY are input data values, refY is reference.
        /// Output files:
        ///    - *SplineData.csv {dataX, dataY} (input data)
        ///    - *SplineValues.csv {dataX, dataY, refY, fitY},
        ///      where fitY is spline value (i.e. approximated value);
        ///    - *SplineCoefficients.csv {dataXleft, dataXright, a, b, c, d},
        ///      where a, b, c, d are spline coefficients.
        /// Note that {dataX, dataY} in input and output files may differ because of 
        /// sorting and merging data points with the same x.</summary>
        [Test]
        public void SplineOptimalSmoothTest()
        {
            ClString name = "SplineOptimalSmoothTest";
            BeginTest(name);

            // Testing data: CallGeomBrownian; Linear, Parabolic, Exponential Data; PayoffButterfly, PayoffCallOption.

            var test = new ApproximationTestsHelper(Log, Context, name, ConsoleOutput);
            foreach (ClString testName in test.GetInputFileNames(GetType().FullName))
            {
                test.StartTest(testName);

                // Read input data.
                ClString fileName = testName + ".csv";
                ClPoint[] data, reference;
                test.ReadCSVInputData(fileName, out data, out reference);

                // Initialise spline and calculate spline coefficients.
                var spline = new ClSplineCubicOptimalSmooth1D(new Spline1DBuilder())
                {
                    OptimSteer =
                    {
                        StartSeed = StartSeed,
                        Method = ClSplineCubicOptimalSmooth1D.OptimizationSteering.UseMethod.CvGoldenSection
                    }
                };
                // TODO Switch to CvGoldenSection, CvScan, NcpScan or NcpGoldenSection optimization.
                spline.LoadData(data);
                test.CalculationTime.Restart();
                spline.Calculate();
                test.CalculationTime.Stop();

                // Store results and complete test.
                test.Print(String.Format("Optimal lambda: {0}", spline.Params.SmoothParam));
                if (PlotResults)
                    test.StoreSpline1D(spline, reference);
                test.Print("Completed " + name + " for " + testName);
            }
            test.FinishTest(PrintPerformance);
        }

        /// <summary>Test for convex spline approximation. (class Cl.ClConvexSpline1D).
        /// Input *.csv files: {dataX, dataY, refY},
        ///       where dataX and dataY are input data values, refY is reference.
        /// Output files:
        ///    - *_SplineData.csv {dataX, dataY} (input data)
        ///    - *SplineValues.csv {dataX, dataY, refY, fitY},
        ///      where fitY is spline value (i.e. approximated value);
        ///    - *_SplineCoefficients.csv {dataXleft, dataXright, a, b, c, d},
        ///      where a, b, c, d are spline coefficients.
        ///    - *_Spline.png (plot)
        /// Note that {dataX, dataY} in input and output files may differ because of 
        /// sorting and merging data points with the same x.</summary>
        [Test]
        public void ConvexSplineTest()
        {
            ClString name = "ConvexSplineTest";
            BeginTest(name);
            const double convexityConstraint = 0.5;
            const double step = 0.1;
            var test = new ConvexSplineTest(Log, Context, name, ConsoleOutput);

            // Testing data: CallGeomBrownian; FunctionExp(X), FunctionXPow2, FunctionXPow4; PayoffCallOption.
            
            foreach (ClString testName in test.GetInputFileNames(GetType().FullName))
            {
                test.StartTest(testName);
                // Read input data.
                ClString fileName = testName + ".csv";
                ClPoint[] data, reference;
                test.ReadCSVInputData(fileName, out data, out reference);
                // Initialise spline and calculate spline coefficients.
                test.Print(String.Format("Calculating spline."));

                var convexSpline = new ClConvexSpline1D
                {
                    ConvexityConstraint = convexityConstraint,
                    Step = step
                };
                convexSpline.LoadData(data);
                convexSpline.Extrapolate = false;
                test.CalculationTime.Restart();
                convexSpline.Calculate();
                test.CalculationTime.Stop();
                var dat = convexSpline.GridData(step);
                var smoothingSpline = new ClSplineCubicSmooth1D(new Spline1DBuilder());
                smoothingSpline.LoadData(data);
                smoothingSpline.Calculate();
                if (PlotResults)
                    test.StoreSpline1D(smoothingSpline, convexSpline, data, reference);
                Log.Result(String.Format(name + " for " + testName + " completed." + "{0}" +
                    "Number of points: " + data.Length + ";{0}" +
                    "Processing time: " + test.CalculationTime.ElapsedMilliseconds + ";{0}" +
                    "Deviation between spline and reference: " + ClApproximUtils.DeviationBetween(convexSpline, reference) + ";{0}" +
                    "Convex measure: " + convexSpline.ConvexMeasure(dat) + ";{0}", Environment.NewLine));
                test.CalculationTime.Reset();
            }
            Log.Result(name + " completed");
            test.FinishTest(PrintPerformance);
        }

        ///<summary>Test for monotone cubic spline implementation. 
        /// (class ClMonotoneSpline1D).</summary>
        [Test]
        public void MonotoneSplineTest()
        {
            ClString name = "MonotoneSplineTest";
            BeginTest(name);
            var test = new MonotoneSplineTest(Log, Context, name, ConsoleOutput);

            // Testing data: CallGeomBrownian; f(x) = x, F(x) = x^4 + 2*x^3 + 5*x*x - 6*x + 7,
            // f(x) = sin(x), f(x) = 1/2 (x-cos(x) sin(x)), f(x) = x+e^-x, etc.
            
            foreach (ClString testName in test.GetInputFileNames(GetType().FullName))
            {
                test.StartTest(testName);
                // Read input data.
                ClString fileName = testName + ".csv";
                ClPoint[] data, reference;
                test.ReadCSVInputData(fileName, out data, out reference);
                // Initialise spline and calculate spline coefficients.
                test.Print(String.Format("Calculating spline."));
                var monotoneSpline = new ClMonotoneSpline1D();
                test.CalculationTime.Restart();
                monotoneSpline.LoadData(data);
                monotoneSpline.Calculate();
                test.CalculationTime.Stop();
                var gridData = monotoneSpline.ValueOnGrid(0.01);
                var dat = new ClPoint[gridData.RowCount];
                for (var i = 0; i < dat.Length; i++)
                    dat[i] = new ClPoint(gridData[i, 0], gridData[i, 1]);
                var smoothingSpline = new ClSplineCubicSmooth1D(new Spline1DBuilder());
                smoothingSpline.LoadData(data);
                smoothingSpline.Calculate();
                if (PlotResults)
                    test.StoreSpline1D(smoothingSpline, monotoneSpline, data, reference);
                Log.Result(String.Format(name + " for " + testName + " completed." + "{0}" +
                    "Number of points: " + data.Length + ";{0}" +
                    "Processing time: " + test.CalculationTime.ElapsedMilliseconds + ";{0}" +
                    "Deviation between spline and reference: " + ClApproximUtils.DeviationBetween(monotoneSpline, reference) + ";{0}", Environment.NewLine));
            }
            Log.Result(name + " completed");
            test.FinishTest(PrintPerformance);
        }

        [Test]
        public void KernelRegressionTest()
        {
            ClString name = "KernelRegressionTest";
            BeginTest(name);

            // Testing data: CallGeomBrownian; Linear, Parabolic, Exponential Data; PayoffButterfly, PayoffCallOption.

            NadarayaWatsonKernelRegressionTest test = new NadarayaWatsonKernelRegressionTest(Log, Context, name, ConsoleOutput);
            foreach (ClString testName in test.GetInputFileNames(GetType().FullName))
            {
                test.StartTest(testName);

                // Read input data.
                ClString fileName = testName + ".csv";
                ClPoint[] data, reference;
                test.ReadCSVInputData(fileName, out data, out reference);
                test.TotalPointsNumber += data.Length;

                // Initialise spline and calculate spline coefficients.
                test.CalculationTime.Restart();
                ClNadarayaWatsonKernelRegression spline = new ClNadarayaWatsonKernelRegression(new ClKernelRegressionParams());
                //spline.OptimSteer.StartSeed = startSeed;
                //// TODO Switch to CvGoldenSection, CvScan, NcpScan or NcpGoldenSection optimization.
                //spline.OptimSteer.Method = ClSplineCubicOptimalSmooth1D.OptimizationSteering.UseMethod.CvGoldenSection;
                spline.LoadData(data);
                //spline.Calculate();
                test.CalculationTime.Stop();
                Array.Sort(data, new ClPointComparerX());
                ClCsvMatrix gridData = test.ValueOnGrid(spline, data[0].X[0], data[data.Length - 1].X[0], 0.01);
                // Store results and complete test.
                if (PlotResults)
                    test.StoreNadarayaWatson1D(spline, data, reference, gridData);
                test.Print("Completed " + name + " for " + testName);
            }
            test.FinishTest(PrintPerformance);
        }

        /// <summary>Test smoothed cubic B-splines with given slopes at end points as boundary conditions (class Cl.ClSplineCubicSmoothSlope1D).
        /// Input:
        ///    - *.csv files: {dataX, dataY, refY},
        ///       where dataX and dataY are input data values, refY is reference.
        ///    - *Slope.csv {left right},
        ///       where left and right are required left and right slopes of the spline, respectively.
        /// Output files:
        ///    - *SplineData.csv {dataX, dataY} (input data)
        ///    - *SplineValues.csv {dataX, dataY, refY, fitY},
        ///      where fitY is spline value (i.e. approximated value);
        ///    - *SplineCoefficients.csv {dataXleft, dataXright, a, b, c, d},
        ///      where a, b, c, d are spline coefficients.
        /// Note that {dataX, dataY} in input and output files may differ because of 
        /// sorting and merging data points with the same x.</summary>
        [Test]
        public void SplineSlopeConditionsTest()
        {
            ClString name = "SplineSlopeConditionsTest";
            BeginTest(name);

            // Testing data: CallGeomBrownian; Linear, Parabolic, Exponential Data; PayoffButterfly, PayoffCallOption.

            ApproximationTestsHelper test = new ApproximationTestsHelper(Log, Context, name, ConsoleOutput);
            foreach (ClString testName in test.GetInputFileNames(GetType().FullName, "Slope"))
            {
                test.StartTest(testName);

                // Read input data.
                ClString fileName = testName + ".csv";
                ClPoint[] data, reference;
                test.ReadCSVInputData(fileName, out data, out reference);

                // Read slope.
                double slopeLeft = 0;
                double slopeRight = 0;
                var n = reference.Length;
                try
                {
                    fileName = testName + "Slope.csv";
                    var matrix = new ClCsvMatrix(Context.AsFileContext.LoadText(fileName));
                    slopeLeft = matrix[0, 0].ToDouble();
                    slopeRight = matrix[0, 1].ToDouble();
                }
                catch (Exception)
                {
                    // if slopes are not set, calculates slopes by first and last two points of reference points.
                    slopeLeft = (reference[1].Value - reference[0].Value) / (reference[1].X[0] - reference[0].X[0]);
                    slopeRight = (reference[n - 1].Value - reference[n - 2].Value) / (reference[n - 1].X[0] - reference[n - 2].X[0]);
                }
                
                // Initialise spline and calculate spline coefficients.
                var spline = new ClSplineCubicSmoothSlope1D(new Spline1DBuilder())
                {
                    Params =
                    {
                        SlopeLeft = slopeLeft,
                        SlopeRight = slopeRight
                    }
                };
                spline.LoadData(data);
                test.CalculationTime.Restart();
                spline.Calculate();
                test.CalculationTime.Stop();
                // Store results and complete test.
                test.Print(String.Format("Optimal lambda: {0}", spline.Params.SmoothParam));
                if (PlotResults)
                    test.StoreSpline1D(spline, reference);
                test.Print("Completed " + name + " for " + testName);
            }
            test.FinishTest(PrintPerformance);
        }

        /// <summary>Test bounded smoothed cubic B-splines with given slopes at end points as boundary conditions (class Cl.ClSplineCubicSmoothSlope1D).
        /// Input:
        ///    - *.csv files: {dataX, dataY, refY},
        ///       where dataX and dataY are input data values, refY is reference.
        ///    - *Slope.csv {left right},
        ///       where left and right are required left and right slopes of the spline, respectively.
        ///    - *Bounds.csv {min max},
        ///       where min and max are required minimum and maximum boundaries of the spline, respectively.
        /// Output files:
        ///    - *SplineData.csv {dataX, dataY} (input data)
        ///    - *SplineValues.csv {dataX, dataY, refY, fitY},
        ///      where fitY is spline value (i.e. approximated value);
        ///    - *SplineCoefficients.csv {dataXleft, dataXright, a, b, c, d},
        ///      where a, b, c, d are spline coefficients.
        /// Note that {dataX, dataY} in input and output files may differ because of 
        /// sorting and merging data points with the same x.</summary>
        [Test]
        public void SplineSlopeBoundTest()
        {
            ClString name = "SplineSlopeBoundTest";
            BeginTest(name);

            // Testing data: CallGeomBrownian; Linear, Parabolic, Exponential Data; PayoffButterfly, PayoffCallOption.

            ApproximationTestsHelper test = new ApproximationTestsHelper(Log, Context, name, ConsoleOutput);
            foreach (ClString testName in test.GetInputFileNames(GetType().FullName, "Slope", "Bounds"))
            {
                test.StartTest(testName);

                // Read input data.
                ClString fileName = testName + ".csv";
                ClPoint[] data, reference;
                test.ReadCSVInputData(fileName, out data, out reference);

                // Read slope.
                // Read slope.
                double slopeLeft = 0;
                double slopeRight = 0;
                var n = reference.Length;
                try
                {
                    fileName = testName + "Slope.csv";
                    var matrix = new ClCsvMatrix(Context.AsFileContext.LoadText(fileName));
                    slopeLeft = matrix[0, 0].ToDouble();
                    slopeRight = matrix[0, 1].ToDouble();
                }
                catch (Exception)
                {
                    // if slopes are not set, calculates slopes by first and last two points of reference points.
                    slopeLeft = (reference[1].Value - reference[0].Value) / (reference[1].X[0] - reference[0].X[0]);
                    slopeRight = (reference[n - 1].Value - reference[n - 2].Value) / (reference[n - 1].X[0] - reference[n - 2].X[0]);
                }

                // Read bounds.
                ClDouble boundMin;
                ClDouble boundMax;
                try
                {
                    fileName = testName + "Bounds.csv";
                    var matrix = new ClCsvMatrix(Context.AsFileContext.LoadText(fileName));
                    boundMin = matrix[0, 0].ToDouble();
                    boundMax = matrix[0, 1].ToDouble();
                }
                catch (Exception)
                {
                    var refer = new double[reference.Length];
                    for (var i = 0; i < refer.Length; i++)
                        refer[i] = reference[i].Value;
                    boundMin = refer.Min();
                    boundMax = refer.Max();
                }
                
                // Initialise spline and calculate spline coefficients.
                var splineParams = new Spline1DBuilder
                {
                    Min = boundMin,
                    Max = boundMax,
                    SlopeLeft = slopeLeft,
                    SlopeRight = slopeRight
                };
                var spline = new ClSplineCubicSmoothSlope1D(splineParams);
                spline.LoadData(data);
                test.CalculationTime.Restart();
                spline.Calculate();
                test.CalculationTime.Stop();


                // Store results and complete test.
                if (PlotResults)
                    test.StoreSpline1D(spline, reference);
                test.Print("Completed " + name + " for " + testName);
            }
            test.FinishTest(PrintPerformance);
        }

        /// <summary>Test smoothed cubic B-splines with fixed values at the end points as boundary conditions (class Cl.ClSplineCubicSmoothEndValue1D).
        /// Input:
        ///    - *.csv files: {dataX, dataY, refY},
        ///       where dataX and dataY are input data values, refY is reference.
        ///    - *Value.csv {left right},
        ///       where left and right are required left and right values of the spline, respectively.
        /// Output files:
        ///    - *SplineData.csv {dataX, dataY} (input data)
        ///    - *SplineValues.csv {dataX, dataY, refY, fitY},
        ///      where fitY is spline value (i.e. approximated value);
        ///    - *SplineCoefficients.csv {dataXleft, dataXright, a, b, c, d},
        ///      where a, b, c, d are spline coefficients.
        /// Note that {dataX, dataY} in input and output files may differ because of 
        /// sorting and merging data points with the same x.</summary>

        [Test]
        public void SplineEndValueTest()
        {
            ClString name = "SplineEndValueTest";
            BeginTest(name);

            // Testing data: Call, Put, Risk Reversal, Straddle, Strangle, Bull, Bear, Butterfly and Barrier (Down-and-in call, Down-and-out call, Up-and-in call, Up-and-out call ) Options.

            ApproximationTestsHelper test = new ApproximationTestsHelper(Log, Context, name, ConsoleOutput);
            foreach (ClString testName in test.GetInputFileNames(GetType().FullName, "Value"))
            {
                test.StartTest(testName);

                // Read input data.
                ClString fileName = testName + ".csv";
                ClPoint[] data, reference;
                test.ReadCSVInputData(fileName, out data, out reference);

                // Read end values.
                double valueLeft = 0;
                double valueRight = 0;
                var n = reference.Length;
                try
                {
                    fileName = testName + "Value.csv";
                    var matrix = new ClCsvMatrix(Context.AsFileContext.LoadText(fileName));
                    valueLeft = matrix[0, 0].ToDouble();
                    valueRight = matrix[0, 1].ToDouble();
                }
                catch (Exception)
                {
                    // set reference values at end points if end values are not set.
                    valueLeft = reference[0].Value;
                    valueRight = reference[n - 1].Value;
                }

                // Initialise spline and calculate spline coefficients.
                Spline1DBuilder builder = new Spline1DBuilder();
                builder.ValueLeft = valueLeft;
                builder.ValueRight = valueRight;
                var spline = new ClSplineCubicSmoothEndValue1D(builder);

                spline.LoadData(data);
                test.CalculationTime.Restart();
                spline.Calculate();
                test.CalculationTime.Stop();

                // Store results and complete test.
                test.Print(String.Format("Optimal lambda: {0}", spline.Params.SmoothParam));
                if (PlotResults)
                    test.StoreSpline1D(spline, reference);
                test.Print("Completed " + name + " for " + testName);
            }
            test.FinishTest(PrintPerformance);
        }

        /// <summary>Test smoothed cubic B-splines with fixed values at the end points as boundary conditions (class Cl.ClSplineCubicSmoothEndValue1D).
        /// Input:
        ///    - *.csv files: {dataX, dataY, refY},
        ///       where dataX and dataY are input data values, refY is reference.
        ///    - *Value.csv {left right},
        ///       where left and right are required left and right values of the spline, respectively.
        ///    - *Slope.csv {left right},
        ///       where left and right are required left and right slopes of the spline, respectively.
        /// Output files:
        ///    - *SplineData.csv {dataX, dataY} (input data)
        ///    - *SplineValues.csv {dataX, dataY, refY, fitY},
        ///      where fitY is spline value (i.e. approximated value);
        ///    - *SplineCoefficients.csv {dataXleft, dataXright, a, b, c, d},
        ///      where a, b, c, d are spline coefficients.
        /// Note that {dataX, dataY} in input and output files may differ because of 
        /// sorting and merging data points with the same x.</summary>

        [Test]
        public void SplineSlopeValueTest()
        {
            ClString name = "SplineSlopeValueTest";
            BeginTest(name);

            // Testing data: CallGeomBrownian; Linear, Parabolic, Exponential Data; PayoffButterfly, PayoffCallOption.

            ApproximationTestsHelper test = new ApproximationTestsHelper(Log, Context, name, ConsoleOutput);
            foreach (ClString testName in test.GetInputFileNames(GetType().FullName, "Slope", "Value"))
            {
                test.StartTest(testName);

                // Read input data.
                ClString fileName = testName + ".csv";
                ClPoint[] data, reference;
                test.ReadCSVInputData(fileName, out data, out reference);

                // Read slope.
                double slopeLeft = 0;
                double slopeRight = 0;
                var n = reference.Length;
                try
                {
                    fileName = testName + "Slope.csv";
                    var matrix = new ClCsvMatrix(Context.AsFileContext.LoadText(fileName));
                    slopeLeft = matrix[0, 0].ToDouble();
                    slopeRight = matrix[0, 1].ToDouble();
                }
                catch (Exception)
                {
                    // if slopes are not set, calculates slopes by first and last two points of reference points.
                    slopeLeft = (reference[1].Value - reference[0].Value) / (reference[1].X[0] - reference[0].X[0]);
                    slopeRight = (reference[n - 1].Value - reference[n - 2].Value) / (reference[n - 1].X[0] - reference[n - 2].X[0]);
                }

                // Read end values.
                double valueLeft = 0;
                double valueRight = 0;
                try
                {
                    fileName = testName + "Value.csv";
                    var matrix = new ClCsvMatrix(Context.AsFileContext.LoadText(fileName));
                    valueLeft = matrix[0, 0].ToDouble();
                    valueRight = matrix[0, 1].ToDouble();
                }
                catch (Exception)
                {
                    // set reference values at end points if end values are not set.
                    valueLeft = reference[0].Value;
                    valueRight = reference[n - 1].Value;
                }

                // Initialise spline and calculate spline coefficients.
                var builder = new Spline1DBuilder
                {
                    ValueLeft = valueLeft,
                    ValueRight = valueRight,
                    SlopeLeft = slopeLeft,
                    SlopeRight = slopeRight
                };
                var spline = new ClSplineCubicSmoothSlopeValue1D(builder);

                spline.LoadData(data);
                test.CalculationTime.Restart();
                spline.Calculate();
                test.CalculationTime.Stop();

                // Store results and complete test.
                test.Print(String.Format("Optimal lambda: {0}", spline.Params.SmoothParam));
                if (PlotResults)
                    test.StoreSpline1D(spline, reference);
                test.Print("Completed " + name + " for " + testName);
            }
            test.FinishTest(PrintPerformance);
        }
    }
}
