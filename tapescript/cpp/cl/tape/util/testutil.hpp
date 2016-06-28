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

#ifndef cl_tape_util_testutil_hpp
#define cl_tape_util_testutil_hpp

#include <vector>
#include <time.h>

// Return vector of random double numbers from range [0, 10).
inline std::vector<double> getRandomVector(int dim)
{
    srand((unsigned int)time(nullptr));
    std::vector<double> result(dim);
    std::for_each(result.begin(), result.end(), [](double & item)
    {
        item = rand() % 10 + (rand() % 1000 / 1000.0);
    });
    return result;
}

// Casts vector elements from RType to LType and returns vector of RType elements.
template <typename LType, typename RType>
inline std::vector<LType> vectorCast(std::vector<RType> const& rhs)
{
    int dim = rhs.size();
    std::vector<LType> result(dim);
    for (int i = 0; i < dim; i++)
    {
        result[i] = (RType)rhs[i];
    }
    return result;
}

inline std::string currentTime()
{
    time_t now = time(nullptr);
    struct tm * tm_today = localtime(&now);

    std::stringstream ss;
    ss << tm_today->tm_hour << "h " << tm_today->tm_min << "m "
        << tm_today->tm_sec << "s ";

    return ss.str();
}

#endif // cl_tape_util_testutil_hpp