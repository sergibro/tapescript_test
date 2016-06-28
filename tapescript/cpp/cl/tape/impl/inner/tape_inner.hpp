/*
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

#ifndef cl_tape_impl_inner_tape_inner_hpp
#define cl_tape_impl_inner_tape_inner_hpp

#include <limits>
#include <valarray>
#include <sstream>
#include <assert.h>

#include <cl/tape/impl/inner/array_traits.hpp>

namespace cl
{
    /// <summary>Class that aggregates a scalar value and an array value.
    /// Used as Base template parameter for CppAD::AD class.</summary>
    template <class Array>
    struct tape_inner
    {
        typedef array_traits<Array> traits;
        typedef typename traits::scalar_type scalar_type;
        typedef typename traits::array_type array_type;
        typedef typename traits::size_type size_type;

        enum Mode
        {
            ScalarMode = 1 << 0
            , IntrusiveScalar = 1 << 1
            , ArrayMode = 1 << 2
        };

        // Default and scalar_type constructor.
        tape_inner(const scalar_type& val = scalar_type())
            : mode_(ScalarMode)
            , scalar_value_(val)
            , array_value_()
        {}

        tape_inner(const tape_inner&) = default;

        tape_inner(tape_inner&& other)
            : mode_(other.mode_)
            , scalar_value_(std::move(other.scalar_value_))
            , array_value_(std::move(other.array_value_))
        {}

        // Array mode is used for array value storage.
        tape_inner(const array_type& v)
            : mode_(ArrayMode)
            , scalar_value_()
            , array_value_(v)
        {}

        tape_inner(array_type&& v)
            : mode_(ArrayMode)
            , scalar_value_()
            , array_value_(std::move(v))
        {}

        // Construct as array with equal coefficients.
        tape_inner(const scalar_type& val, size_t count)
            : mode_(ArrayMode)
            , scalar_value_()
            , array_value_(traits::make(val, count))
        {}

        // Construct as array with values passed by pointer.
        tape_inner(const scalar_type* ptr, size_t count)
            : mode_(ArrayMode)
            , scalar_value_()
            , array_value_(traits::make(ptr, count))
        {}

        // Construct as array with values passed by initializer_list.
        tape_inner(std::initializer_list<scalar_type> il)
            : tape_inner(il.begin(), il.size())
        {}

        inline tape_inner& operator=(tape_inner&& other)
        {
            mode_ = other.mode_;
            scalar_value_ = std::move(other.scalar_value_);
            array_value_ = std::move(other.array_value_);
            return *this;
        }

        // Returns true if scalar mode used (ordinary or intrusive).
        bool is_scalar() const
        {
            return (mode_ & (ScalarMode | IntrusiveScalar)) != 0;
        }

        // Returns true if array mode used.
        bool is_array() const
        {
            return !is_scalar();
        }

        // Returns true if intrusive scalar mode used.
        bool is_intrusive() const
        {
            return mode_ == IntrusiveScalar;
        }

        void set_intrusive()
        {
            if (mode_ == ScalarMode)
            {
                mode_ = IntrusiveScalar;
            }
        }

        void set_not_intrusive()
        {
            if (mode_ == IntrusiveScalar)
            {
                mode_ = ScalarMode;
            }
        }

        // Converts to scalar value.
        scalar_type to_scalar() const
        {
            if (is_array())
            {
                cl::throw_("Not a scalar");
            }
            return scalar_value_;
        }

        // Returns arithmetic negation.
        inline tape_inner operator-() const
        {
            if (is_scalar())
            {
                return -scalar_value_;
            }
            return tape_inner(-array_value_);
        }

        // Returns a new tape_inner of the same size with values which are acquired
        // by applying function func to the previous values of the elements.
        tape_inner apply(scalar_type func(scalar_type)) const
        {
            if (is_scalar())
            {
                return func(scalar_value_);
            }
            array_type result(array_value_.size());
            for (size_type i = 0; i < result.size(); i++)
            {
                result[i] = func(array_value_[i]);
            }
            return result;
        }

        // Assign operations.
#define CL_INNER_ARRAY_ASSIGN_OPERATOR(Op)                                                      \
        inline tape_inner& operator Op ## = (const tape_inner& right)                           \
        {                                                                                       \
            if (is_intrusive())                                                                 \
            {                                                                                   \
                scalar_value_ Op##= right.sum();                                                \
            }                                                                                   \
            else if (is_scalar() && right.is_scalar())                                          \
            {                                                                                   \
                scalar_value_ Op##= right.scalar_value_;                                        \
            }                                                                                   \
            else if (is_array() && right.is_scalar())                                           \
            {                                                                                   \
                array_value_ Op##= right.scalar_value_;                                         \
            }                                                                                   \
            else if (is_scalar() && right.is_array())                                           \
            {                                                                                   \
                array_value_ = scalar_value_ Op right.array_value_;                             \
                mode_ = ArrayMode;                                                              \
            }                                                                                   \
            else if (is_array() && right.is_array())                                            \
            {                                                                                   \
                array_value_ Op##= right.array_value_;                                          \
            }                                                                                   \
            return *this;                                                                       \
        }
        CL_INNER_ARRAY_ASSIGN_OPERATOR(+)
        CL_INNER_ARRAY_ASSIGN_OPERATOR(-)
        CL_INNER_ARRAY_ASSIGN_OPERATOR(*)
        CL_INNER_ARRAY_ASSIGN_OPERATOR(/)
#undef CL_INNER_ARRAY_ASSIGN_OPERATOR

        // Gives the element of an array by the index.
        // If the object is array valued returns an element of the array value.
        // Othervise returns scalar value.
        // Use this method for element-wise operations.
        const scalar_type& element_at(size_t index) const
        {
            if (is_scalar())
            {
                return scalar_value_;
            }
            return array_value_[(size_type)index];
        }

        scalar_type sum() const
        {
            if (is_scalar())
            {
                return scalar_value_;
            }
            scalar_type result = 0.0;
            for (size_t i = 0; i < size(); i++)
            {
                result += array_value_[i];
            }
            return result;
        }

        // Returns size of array value.
        size_t size() const
        {
            CL_ASSERT(is_array(), "Have to be an array.");
            size_type temp = array_value_.size();
            CL_ASSERT(temp >= 0, "");
            return (size_t)temp;
        }

        void resize(size_t size)
        {
            mode_ = ArrayMode;
            array_value_.resize(size);
        }

        scalar_type* begin()
        {
            CL_ASSERT(is_array(), "Have to be an array to take an iterator.");
            return &array_value_[0];
        }

        scalar_type const* begin() const
        {
            CL_ASSERT(is_array(), "Have to be an array to take an iterator.");
            return &array_value_[0];
        }

        scalar_type* end()
        {
            CL_ASSERT(is_array(), "Have to be an array to take an iterator.");
            return &array_value_[0] + size();
        }

        scalar_type const* end() const
        {
            CL_ASSERT(is_array(), "Have to be an array to take an iterator.");
            return &array_value_[0] + size();
        }

        scalar_type& operator[](size_t index)
        {
            CL_ASSERT(is_array(), "Have to be an array to access via [].");
            return array_value_[index];
        }

        scalar_type const& operator[](size_t index) const
        {
            CL_ASSERT(is_array(), "Have to be an array to access via [].");
            return array_value_[index];
        }

        Mode mode_;
        scalar_type scalar_value_;
        array_type array_value_;
    };


    // Stream insertion operator.
    template <class Array>
    inline std::ostream& operator<<(std::ostream& os, const tape_inner<Array>& x)
    {
        if (x.is_scalar())
        {
            return os << x.scalar_value_;
        }
        if (x.size() == 0)
        {
            return os << "{}";
        }

        std::stringstream ss;
        ss.precision(os.precision());
        ss << "{ " << x.array_value_[0];
        for (size_t i = 1; i < x.size(); ++i)
        {
            ss << ", " << x.array_value_[i];
        }
        ss << " }";

        return os << ss.str();
    }


    // Arithmetic binary operations.
#define CL_BIN_INNER_ARRAY_OPERATOR(Op)                                                         \
    template <class Array>                                                                      \
    inline tape_inner<Array> operator Op(                                                       \
        const tape_inner<Array>& x                                                              \
        , const tape_inner<Array>& y)                                                           \
    {                                                                                           \
        if (x.is_scalar() && y.is_scalar())                                                     \
        {                                                                                       \
            return x.scalar_value_ Op y.scalar_value_;                                          \
        }                                                                                       \
        else if (x.is_array() && y.is_scalar())                                                 \
        {                                                                                       \
            return tape_inner<Array>(x.array_value_ Op y.scalar_value_);                        \
        }                                                                                       \
        else if (x.is_scalar() && y.is_array())                                                 \
        {                                                                                       \
            return tape_inner<Array>(x.scalar_value_ Op y.array_value_);                        \
        }                                                                                       \
        else /* (x.is_array() && y.is_array()) */                                               \
        {                                                                                       \
            return tape_inner<Array>(x.array_value_ Op y.array_value_);                         \
        }                                                                                       \
    }                                                                                           \
                                                                                                \
    template <class Array>                                                                      \
    inline tape_inner<Array> operator Op(                                                       \
        const tape_inner<Array>& x                                                              \
        , const typename tape_inner<Array>::scalar_type& y)                                     \
    {                                                                                           \
        if (x.is_scalar())                                                                      \
        {                                                                                       \
            return x.scalar_value_ Op y;                                                        \
        }                                                                                       \
        return tape_inner<Array>(x.array_value_ Op y);                                          \
    }                                                                                           \
                                                                                                \
    template <class Array>                                                                      \
    inline tape_inner<Array> operator Op(                                                       \
        const typename tape_inner<Array>::scalar_type& x                                        \
        , const tape_inner<Array>& y)                                                           \
    {                                                                                           \
        if (y.is_scalar())                                                                      \
        {                                                                                       \
            return x Op y.scalar_value_;                                                        \
        }                                                                                       \
        return tape_inner<Array>(x Op y.array_value_);                                          \
    }

    CL_BIN_INNER_ARRAY_OPERATOR(-)
    CL_BIN_INNER_ARRAY_OPERATOR(*)
    CL_BIN_INNER_ARRAY_OPERATOR(/)
    CL_BIN_INNER_ARRAY_OPERATOR(+)
#undef CL_BIN_INNER_ARRAY_OPERATOR


    // Logical binary operations.
#define CL_BOOL_INNER_ARRAY_OPERATOR(Op, Code)                                                  \
    template <class Array>                                                                      \
    inline bool operator Op(const tape_inner<Array>& x, const tape_inner<Array>& y)             \
    {                                                                                           \
        if (x.is_scalar() && y.is_scalar())                                                     \
        {                                                                                       \
            return x.scalar_value_ Op y.scalar_value_;                                          \
        }                                                                                       \
        bool result = true;                                                                     \
        if (x.is_array() && y.is_scalar())                                                      \
        {                                                                                       \
            result = array_traits<Array>::operator_##Code(x.array_value_, y.scalar_value_);     \
        }                                                                                       \
        else if (x.is_scalar() && y.is_array())                                                 \
        {                                                                                       \
            result = array_traits<Array>::operator_##Code(x.scalar_value_, y.array_value_);     \
        }                                                                                       \
        else /* (left.is_array() && right.is_array()) */                                        \
        {                                                                                       \
            result = array_traits<Array>::operator_##Code(x.array_value_, y.array_value_);      \
        }                                                                                       \
        return result;                                                                          \
    }                                                                                           \
                                                                                                \
    template <class Array>                                                                      \
    inline bool operator Op(                                                                    \
        const tape_inner<Array>& x                                                              \
        , const typename tape_inner<Array>::scalar_type& y)                                     \
    {                                                                                           \
        if (x.is_scalar())                                                                      \
        {                                                                                       \
            return x.scalar_value_ Op y;                                                        \
        }                                                                                       \
        return array_traits<Array>::operator_##Code(x.array_value_, y);                         \
    }                                                                                           \
                                                                                                \
    template <class Array>                                                                      \
    inline bool operator Op(                                                                    \
        const typename tape_inner<Array>::scalar_type& x                                        \
        , const tape_inner<Array>& y)                                                           \
    {                                                                                           \
        if (y.is_scalar())                                                                      \
        {                                                                                       \
            return x Op y.scalar_value_;                                                        \
        }                                                                                       \
        return array_traits<Array>::operator_##Code(x, y.array_value_);                         \
    }

    CL_BOOL_INNER_ARRAY_OPERATOR(!=, Ne)
    CL_BOOL_INNER_ARRAY_OPERATOR(==, Eq)
    CL_BOOL_INNER_ARRAY_OPERATOR(< , Lt)
    CL_BOOL_INNER_ARRAY_OPERATOR(<=, Le)

    template <class Array>
    inline bool operator>(const tape_inner<Array>& x, const tape_inner<Array>& y)
    {
        return y < x;
    }
    template <class Array>
    inline bool operator>(const tape_inner<Array>& x, const typename tape_inner<Array>::scalar_type& y)
    {
        return y < x;
    }
    template <class Array>
    inline bool operator>(const typename tape_inner<Array>::scalar_type& x, const tape_inner<Array>& y)
    {
        return y < x;
    }
    template <class Array>
    inline bool operator>=(const tape_inner<Array>& x, const tape_inner<Array>& y)
    {
        return y <= x;
    }
    template <class Array>
    inline bool operator>=(const tape_inner<Array>& x, const typename tape_inner<Array>::scalar_type& y)
    {
        return y <= x;
    }
    template <class Array>
    inline bool operator>=(const typename tape_inner<Array>::scalar_type& x, const tape_inner<Array>& y)
    {
        return y <= x;
    }
#undef CL_BOOL_INNER_ARRAY_OPERATOR


    namespace tapescript
    {
        // Standart math functions.
#define CL_INNER_ARRAY_FUNCTION(Name)                                                           \
        template <class Array>                                                                  \
        inline cl::tape_inner<Array> Name(const cl::tape_inner<Array>& x)                       \
        {                                                                                       \
            if (x.is_scalar())                                                                  \
            {                                                                                   \
                return std::Name(x.scalar_value_);                                              \
            }                                                                                   \
            return cl::tape_inner<Array>::traits::Name(x.array_value_);                         \
        }
        CL_INNER_ARRAY_FUNCTION(abs)
        CL_INNER_ARRAY_FUNCTION(acos)
        CL_INNER_ARRAY_FUNCTION(sqrt)
        CL_INNER_ARRAY_FUNCTION(asin)
        CL_INNER_ARRAY_FUNCTION(atan)
        CL_INNER_ARRAY_FUNCTION(cos)
        CL_INNER_ARRAY_FUNCTION(sin)
        CL_INNER_ARRAY_FUNCTION(cosh)
        CL_INNER_ARRAY_FUNCTION(sinh)
        CL_INNER_ARRAY_FUNCTION(exp)
        CL_INNER_ARRAY_FUNCTION(log)
        CL_INNER_ARRAY_FUNCTION(tan)
        CL_INNER_ARRAY_FUNCTION(tanh)
#undef CL_INNER_ARRAY_FUNCTION

        template <class Array>
        inline cl::tape_inner<Array> sign(const cl::tape_inner<Array>& x)
        {
            typedef typename cl::tape_inner<Array>::scalar_type scalar_type;

            if (x.is_scalar())
            {
                if (x.scalar_value_ > 0.)
                    return scalar_type(1.0);
                if (x.scalar_value_ == 0.)
                    return scalar_type(0.0);
                return scalar_type(-1.0);
            }

            return cl::tape_inner<Array>::traits::sign(x.array_value_);
        }

        // Math power functioon.
        template <class Array>
        inline cl::tape_inner<Array> pow(
            const cl::tape_inner<Array>& left
            , const cl::tape_inner<Array>& right)
        {
            typedef typename cl::tape_inner<Array>::traits traits;
            if (left.is_scalar() && right.is_scalar())
            {
                return std::pow(left.scalar_value_, right.scalar_value_);
            }
            else if (left.is_array() && right.is_scalar())
            {
                return traits::pow(left.array_value_, right.scalar_value_);
            }
            else if (left.is_scalar() && right.is_array())
            {
                return traits::pow(left.scalar_value_, right.array_value_);
            }
            else // (left.is_array() && right.is_array())
            {
                return traits::pow(left.array_value_, right.array_value_);
            }
        }

        // Math power functioon.
        template <class Array>
        inline cl::tape_inner<Array> pow(
            const cl::tape_inner<Array>& left
            , const typename cl::tape_inner<Array>::scalar_type& right)
        {
            typedef typename cl::tape_inner<Array>::traits traits;
            if (left.is_scalar())
            {
                return std::pow(left.scalar_value_, right);
            }
            return traits::pow(left.array_value_, right);
        }

        // Math power functioon.
        template <class Array>
        inline cl::tape_inner<Array> pow(
            const typename cl::tape_inner<Array>::scalar_type& left
            , const cl::tape_inner<Array>& right)
        {
            typedef typename cl::tape_inner<Array>::traits traits;
            if (right.is_scalar())
            {
                return std::pow(left, right.scalar_value_);
            }
            return traits::pow(left, right.array_value_);
        }

        template <class T>
        void set_intrusive(T& val, const T& model = T()){}

        template <class Array>
        void set_intrusive(tape_inner<Array>& val, const tape_inner<Array>& model = tape_inner<Array>())
        {
            if (model.is_scalar())
            {
                val.set_intrusive();
            }
        }

        template <class T>
        void set_not_intrusive(T& val){}

        template <class Array>
        void set_not_intrusive(tape_inner<Array>& val)
        {
            val.set_not_intrusive();
        }
    } // namespace tapescript
} // namespace cl


namespace std
{
    template <class Array>
    cl::tape_inner<Array> ceil(cl::tape_inner<Array> const& x)
    {
        return x.apply(std::ceil);
    }

    template <class Array>
    cl::tape_inner<Array> floor(cl::tape_inner<Array> const& x)
    {
        return x.apply(std::floor);
    }

    // CLASS numeric_limits<cl::tape_inner<Array>>
    template<class Array>
    class numeric_limits<cl::tape_inner<Array>>
    {
        typedef typename cl::tape_inner<Array>::scalar_type scalar_type;
    public:
        typedef cl::tape_inner<Array> _Ty;

        static _Ty min()
        {    // return minimum value
            return numeric_limits<scalar_type>::min();
        }

        static _Ty max()
        {    // return maximum value
            return numeric_limits<scalar_type>::max();
        }

        static _Ty epsilon()
        {    // return smallest effective increment from 1.0
            return numeric_limits<scalar_type>::epsilon();
        }
    };
}

#include <cl/tape/impl/inner/tape_inner_fields.hpp>

#endif // cl_tape_impl_inner_tape_inner_hpp
