﻿/*
Copyright (C) 2015-present CompatibL

Performance test results and finance-specific examples are available at:

http://www.tapescript.org

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

#ifndef cl_tape_impl_doublelimits_hpp
#define cl_tape_impl_doublelimits_hpp

#include <cl/tape/impl/double.hpp>
#include <limits>

namespace std
{
    /// <summary>Provides numeric limits for tape_double.</summary>
    template<class Inner>
    class numeric_limits<cl::tape_wrapper<Inner>>
    {
    public:
        static const bool is_specialized = true;

        static Inner min() throw() { return numeric_limits<double>::min(); }
        static Inner max() throw() { return numeric_limits<double>::max(); }
        static Inner lowest() throw() { return numeric_limits<double>::lowest(); }

        static const int digits = numeric_limits<double>::digits;
        static const int digits10 = numeric_limits<double>::digits10;
        static const int max_digits10 = numeric_limits<double>::max_digits10;

        static const bool is_signed = numeric_limits<double>::is_signed;
        static const bool is_integer = numeric_limits<double>::is_integer;
        static const bool is_exact = numeric_limits<double>::is_exact;
        static const int radix = numeric_limits<double>::radix;
        static const Inner epsilon() throw() { return numeric_limits<double>::epsilon(); }
        static const Inner round_error() throw() { return numeric_limits<double>::round_error(); }

        static const int min_exponent = numeric_limits<double>::min_exponent;
        static const int min_exponent10 = numeric_limits<double>::min_exponent10;
        static const int max_exponent = numeric_limits<double>::max_exponent;
        static const int max_exponent10 = numeric_limits<double>::max_exponent10;

        static const bool has_infinity = numeric_limits<double>::has_infinity;
        static const bool has_quiet_NaN = numeric_limits<double>::has_quiet_NaN;
        static const bool has_signaling_NaN = numeric_limits<double>::has_signaling_NaN;

        static Inner infinity() throw() { return numeric_limits<double>::infinity(); }
        static Inner quiet_NaN() throw() { return numeric_limits<double>::quiet_NaN(); }
        static Inner signaling_NaN() throw() { return numeric_limits<double>::signaling_NaN(); }
        static Inner denorm_min() throw() { return numeric_limits<double>::denorm_min(); }
    };
}

#endif  // cl_tape_impl_doublelimits_hpp
