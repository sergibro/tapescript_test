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
    /// <summary> Cubic smoothing spline with slopes and values at the end points as boundary conditions.</summary>
    public class ClSplineCubicSmoothSlopeValue1D : ClSplineCubicSmooth1D
    {
        public ClSplineCubicSmoothSlopeValue1D(Spline1DBuilder parameters) : base(parameters)
        {
            splineMath_ = new SplineCubicSmoothSlopeValue1DUtils();
        }
            
        /// <summary>Calculates spline parameters.</summary>
        public override void Calculate()
        {
            CalculateWithoutExtension();
        }
    }
}
