/*
Copyright (C) 2003-2015 CompatibL

This file is part of TapeScript, an open source library and tape encoding
standard for adjoint algorithmic differentiation (AAD), available from

http://github.com/compatibl/tapescript (source)
http://tapescript.org (documentation)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef cl_tape_examples_impl_spline_cubic_smooth_end_value_hpp
#define cl_tape_examples_impl_spline_cubic_smooth_end_value_hpp

#include <cl/tape/tape.hpp>
#include "impl/utils.hpp"
#include "impl/smooth_cubic_spline.hpp"
#include "impl/spline_cubic_smooth_end_value_utils.hpp"

namespace cl
{
	/// <summary> Cubic smoothing spline with values at the end points as boundary conditions.</summary>
	class spline_cubic_smooth_end_value : public smooth_cubic_spline
	{
	public:
		spline_cubic_smooth_end_value(tobject valueLeft, tobject valueRight)
		{
			spline_math = new spline_cubic_smooth_end_value_utils(valueLeft, valueRight);
		}
	};
}
#endif // cl_tape_examples_impl_spline_cubic_smooth_end_value_hpp