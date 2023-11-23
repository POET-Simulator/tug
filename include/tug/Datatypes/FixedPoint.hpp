#ifndef TUG_FIXED_TYPE
#define TUG_FIXED_TYPE

#include <cstdint>
#include <ostream>
#include <type_traits>

namespace tug {

template <typename BaseType, typename GreaterType, unsigned int FractionBits,
          unsigned int IntegerMult = BaseType(1) << FractionBits>
class FixedPoint {
public:
  // zero init
  constexpr FixedPoint() : _internal_value(0) {}

  // construct a fixed point from value
  template <typename T>
  constexpr FixedPoint(T val)
      : _internal_value(static_cast<BaseType>(val * IntegerMult)) {}

  // cast a fixed point to a data type
  template <typename T> constexpr operator T() const {
    return static_cast<T>(static_cast<T>(_internal_value) / IntegerMult);
  }

  // addition
  template <typename T> constexpr FixedPoint operator+(const T &rhs) const {
    return FixedPoint(add(this->_internal_value, rhs), {});
  }

  // subtraction
  template <typename T> constexpr FixedPoint operator-(const T &rhs) const {
    return FixedPoint(sub(this->_internal_value, rhs), {});
  }

  // multiplication
  template <typename T> constexpr FixedPoint operator*(const T &rhs) const {
    return FixedPoint(mul(this->_internal_value, rhs), {});
  }

  // division
  template <typename T> constexpr FixedPoint operator/(const T &rhs) const {
    return FixedPoint(div(this->_internal_value, rhs), {});
  }

  // equality
  template <typename T> constexpr bool operator==(const T &rhs) const {
    return eq(this->_internal_value, rhs);
  }

private:
  struct placeholder {};
  constexpr FixedPoint(BaseType _val, placeholder) : _internal_value(_val) {}

  static constexpr BaseType add(const BaseType _lhs, const FixedPoint _rhs) {
    return static_cast<BaseType>(_lhs + _rhs._internal_value);
  }

  static constexpr BaseType sub(const BaseType _lhs, const FixedPoint _rhs) {
    return static_cast<BaseType>(_lhs - _rhs._internal_value);
  }

  static constexpr BaseType mul(const BaseType _lhs, const FixedPoint _rhs) {
    return static_cast<BaseType>(
        static_cast<GreaterType>(_lhs * _rhs._internal_value) / IntegerMult);
  }

  static constexpr BaseType div(const BaseType _lhs, const FixedPoint _rhs) {
    return static_cast<BaseType>(static_cast<GreaterType>(_lhs * IntegerMult) /
                                 _rhs._internal_value);
  }

  static constexpr bool eq(const BaseType _lhs, const FixedPoint _rhs) {
    return _lhs == _rhs._internal_value;
  }

  BaseType _internal_value;
};

template <typename B, typename G, unsigned int F>
std::ostream &operator<<(std::ostream &os, const FixedPoint<B, G, F> &fp) {
  os << static_cast<double>(fp);

  return os;
}

template <unsigned int F>
using Fixed32 = FixedPoint<std::int32_t, std::int64_t, F>;

template <unsigned int F>
using Fixed16 = FixedPoint<std::int16_t, std::int32_t, F>;

template <unsigned int F>
using Fixed8 = FixedPoint<std::int8_t, std::int16_t, F>;
} // namespace tug

#endif // TUG_FIXED_TYPE