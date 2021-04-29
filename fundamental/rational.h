/*{ "clmagic/fundamental/rational":{
  "Author":"LongJiangnan",
  "Date":"2019-2021",
  "License":"Free"
} }*/
#pragma once

#include <numeric>
#include <math.h>
namespace calculation {

	template<size_t _Exp, typename _Ti>
	_Ti _Gcdx_approx(_Ti _Ax, _Ti _Bx) {// computes greatest common divisor of _Ax and _Bx
		if ((_Ax >> _Exp) >= _Bx) {
			return (_Ax);
		} else {
			return _Gcdx_approx<_Exp>(_Bx, _Ax % _Bx);
		}
		/*
		Ax     Bx
		16667 10000
		10000 6667
		6667  3333
		3333  1

		16666 10000
		10000 6666
		6666  3334
		3334  3332
		3332  2
		This is simple test
		Result is _Ax when (_Ax modulus _Bx) is very small
		But _Bx must be large than (_Ax mudulus _Bx)
		And not error, because final-divides is <trunc (_Ax modulus _Bx)>
		*/
	}

	template<size_t _Exp, typename _Ti> inline 
	_Ti gcd_approx(_Ti _Ax, _Ti _Bx) {
		if (_Ax == 0 && _Bx == 0) {
			return (1);
		}

		if (_Bx > _Ax) {
			std::swap(_Ax, _Bx);
		}
		return _Gcdx_approx<_Exp>(_Ax < 0 ? -_Ax : _Ax, _Bx < 0 ? -_Bx : _Bx);
	}
	
	template<typename T>
	class rational_ {
		static_assert(std::_Is_any_of_v<T, int32_t, int64_t, intmax_t>, "rational_<integer>");
	public:
		using integer = T;

		constexpr rational_() : numerator(0), denominator(0) {}
		
		constexpr rational_(integer _Nx, integer _Dx = 1) : numerator(_Nx), denominator(_Dx) {}
		
		rational_& operator=(const rational_&) = default;
		
		rational_& operator=(integer _Nx) {
			numerator = _Nx;
			denominator = 1; 
			return *this;
		}
		
		bool operator==(rational_ right) const {
			return this->numerator * right.denominator 
				== right.numerator * this->denominator;
		}
		
		bool operator!=(rational_ right) const {
			return !((*this) == right);
		}
		
		bool operator<(rational_ right) const {
			return this->numerator * right.denominator 
				< right.numerator * this->denominator;
		}
		
		bool operator>(rational_ right) const {
			return right < (*this);
		}
		
		bool operator>=(rational_ right) const {
			return !(*this < right);
		}
		
		bool operator<=(rational_ right) const {
			return !(*this > right);
		}
		
		rational_ operator-() const {
			return rational_(-numerator, denominator);
		}
		
		rational_ operator+(rational_ right) const {
			if ( this->numerator == 0 ) {
				return right;
			}
			if ( right.numerator == 0 ){
				return *this;
			}
			
			/* {
			N1   N2   N1*D2 + N2*D1
			-- + -- = -------------
			D1   D2      D1*D2
			} */
			integer common_divisor = std::gcd(this->denominator, right.denominator);
			integer this_denominator = this->denominator / common_divisor;
			integer right_denominator = right.denominator / common_divisor;
			return rational_(
				this->numerator * right_denominator + right.numerator * this_denominator,
				this_denominator * right_denominator
				);
		}
		
		rational_ operator-(rational_ right) const {
			if ( this->numerator == 0 ) {
				return -right;
			}
			if ( right.numerator == 0 ){
				return *this;
			}

			/* {
			N1   N2   N1*D2 - N2*D1
			-- - -- = -------------
			D1   D2      D1*D2
			} */
			integer common_divisor = std::gcd(this->denominator, right.denominator);
			integer this_denominator = this->denominator / common_divisor;
			integer right_denominator = right.denominator / common_divisor;
			return rational_(
				this->numerator * right_denominator - right.numerator * this_denominator,
				this_denominator * right_denominator
				);
		}
		
		rational_ operator*(rational_ right) const {
			if (this->numerator == 0 || right.numerator == 0) {
				return rational_(0, 0);
			}
			
			/* {
			N1   N2   N1*N2
			-- * -- = -----
			D1   D2   D1*D2
			} */
			integer common_divisor1 = std::gcd(this->numerator, right.denominator);
			integer this_numerator = this->numerator / common_divisor1;
			integer right_denominator = right.denominator / common_divisor1;

			integer common_divisor2 = std::gcd(right.numerator, this->denominator);
			integer right_numerator = right.numerator / common_divisor2;
			integer this_denominator = this->denominator / common_divisor2;
			
			return rational_(
				this_numerator * right_numerator,
				this_denominator * right_denominator
				);
		}

		rational_ operator/(rational_ right) const {
			return (*this) * rational_(right.denominator, right.numerator);
		}
	
		rational_ operator+(integer value) const {
			return rational_(
				numerator + value * denominator, 
				denominator
				);
		}

		rational_ operator-(integer value) const {
			return rational_(
				numerator - value * denominator, 
				denominator
				);
		}

		rational_ operator*(integer value) const {
			integer common_divisor = std::gcd(value, denominator);
			return rational_(
				numerator * (value / common_divisor),
				denominator / common_divisor
				);
		}

		rational_ operator/(integer value) const {
			return (*this) * rational_(1, value);
		}

		rational_& operator++() {
			numerator += denominator;
			return *this;
		}
		
		rational_ operator++(int) {
			rational_ copied = *this;
			++(*this);
			return copied;
		}
		
		rational_& operator+=(rational_ right) {
			return *this = *this + right;
		}
		
		rational_& operator-=(rational_ right) {
			return *this = *this - right;
		}
		
		rational_& operator*=(rational_ right) {
			return *this = *this * right;
		}
		
		rational_& operator/=(rational_ right) {
			return *this = *this / right;
		}
				
		rational_& operator+=(integer value) {
			return (*this) = (*this) + value;
		}
	
		rational_& operator-=(integer value) {
			return (*this) = (*this) - value;
		}
		
		rational_& operator*=(integer value) {
			return (*this) = (*this) * value;
		}
		
		rational_& operator/=(integer value) {
			return (*this) = (*this) / value;
		}
		
		rational_& operator%=(integer value) {
			return (*this) = (*this) % value;
		}

		void reduce() {
			integer common_divisor = std::gcd(numerator, denominator);
			numerator /= common_divisor;
			denominator /= common_divisor;
		}

		/* <formula>
			<source>0 00000001 10101101011010110101000</source>
			<first>
				<tips> IEEE754-floating-formula: (-1)^S * (1+0.Fraction) * 2^(Exponent-Bias) </tips>
					
				(-1)^0 * (1 + 0.10101101011010110101000) * 2^(00000001 - Bias)
				= 1 * 1.10101101011010110101000 * pow(2, _Exp)
				= 1 * 0.110101101011010110101000 * pow(2, _Exp)
				= 1 * 110101101011010110101000/pow(2,_Mn) * pow(2, _Exp)
			</first>				
			<second>
				<tips> pow(2, X) = (1 << X) </tips>

				_Nx     110101101011010110101000
				---- = -------------------------- * ( 1 << _Exp )
				_Dx           1 << _Mn

					110101101011010110101000 << 1
				----------------------------------- * ( 1 << (_Exp - 1) )
							1 << _Mn

					110101101011010110101000
				----------------------------------- * ( 1 << (_Exp - 1) )
							1 << (_Mn-1)

					110101101011010110101000
				----------------------------------- * (1 << 0)
							1 << (_Mn - _Exp)
			</second>
		</formula> */
	
		explicit 
		rational_(float value) {
			constexpr unsigned int sign_mask 
				= 0b10000000000000000000000000000000;
			constexpr unsigned int exponent_mask
				= 0b01111111100000000000000000000000;
			constexpr unsigned int mantissa_mask
				= 0b00000000011111111111111111111111;
			constexpr unsigned int hidden_significant
				= 0b00000000100000000000000000000000;
			constexpr char exp2_bias = 127;

			unsigned int value_bits = reinterpret_cast<uint32_t&>(value);
			unsigned int exp2_bits = (value_bits & exponent_mask) >> 23;
			unsigned int signifi_bits = value_bits & mantissa_mask | hidden_significant;
			char exp2 = reinterpret_cast<char&>(exp2_bits) - exp2_bias;
			exp2 -= 23;
			
			// *this = significant * 2^23
			this->numerator = static_cast<integer>(signifi_bits);
			this->denominator = static_cast<integer>(1);

			// *this = significant * pow(2,exp2)
			if (exp2 > 0) {
				this->numerator = static_cast<integer>(signifi_bits << exp2);
				this->denominator = static_cast<integer>(1);
			} else if (exp2 < 0) {
				static_assert(sizeof(integer) >= sizeof(unsigned int), "ratianal(float)");
				this->numerator = static_cast<integer>(signifi_bits);
				this->denominator = static_cast<integer>(1) << (-exp2);
			} else {
				this->numerator = static_cast<integer>(signifi_bits);
				this->denominator = static_cast<integer>(1);
			}

			// *this *= ~sign
			if ((value_bits & sign_mask) != uint32_t(0)) {
				this->numerator = -this->numerator;
			}

			// divide greater_common_divisor
			integer common_divisor = std::gcd(numerator, denominator);
			numerator /= common_divisor;
			denominator /= common_divisor;
		}

		explicit 
		rational_(double value) {
			constexpr unsigned long long sign_mask     
				= 0b1000000000000000000000000000000000000000000000000000000000000000;
			constexpr unsigned long long exponent_mask
				= 0b0111111111110000000000000000000000000000000000000000000000000000;
			constexpr unsigned long long mantissa_mask
				= 0b0000000000001111111111111111111111111111111111111111111111111111;
			constexpr unsigned long long hidden_significant
				= 0b0000000000010000000000000000000000000000000000000000000000000000;
			constexpr short exp2_bias = 1023;

			// seperate bits
			unsigned long long value_bits = reinterpret_cast<unsigned long long&>(value);
			unsigned long long exp2_bits = (value_bits & exponent_mask) >> 52;
			unsigned long long signifi_bits = value_bits & mantissa_mask | hidden_significant;
			short exp2 = reinterpret_cast<short&>(exp2_bits) - exp2_bias;
			exp2 -= 52;
			
			// *this = significant * pow(2,exp2)
			if (exp2 > 0) {
				this->numerator = static_cast<integer>(signifi_bits << exp2);
				this->denominator = static_cast<integer>(1);
			} else if (exp2 < 0) {
				static_assert(sizeof(integer) >= sizeof(unsigned long long), "ratianal(double)");
				this->numerator = static_cast<integer>(signifi_bits);
				this->denominator = static_cast<integer>(1) << (-exp2);
			} else {
				this->numerator = static_cast<integer>(signifi_bits);
				this->denominator = static_cast<integer>(1);
			}

			// *this *= ~sign
			if ((value_bits & sign_mask) != 0) {
				this->numerator = -this->numerator;
			}

			// divide greater_common_divisor
			integer common_divisor = std::gcd(numerator, denominator);
			numerator /= common_divisor;
			denominator /= common_divisor;
		}

		integer numerator;
		integer denominator;
	};

	using rational32 = rational_<int>;
	using rational64 = rational_<long long>;

	template<typename T> inline
	rational_<T> sin(rational_<T> x) {
		if constexpr (sizeof(T) <= sizeof(float)) {
			return rational_<T>(
				_CSTD sinf(static_cast<float>(x.numerator) / static_cast<float>(x.denominator))
				);
		} else {
			return rational_<T>(
				_CSTD sin(static_cast<double>(x.numerator) / static_cast<double>(x.denominator))
				);
		}
	}

	template<typename T> inline
	rational_<T> cos(rational_<T> x) {
		if constexpr (sizeof(T) <= sizeof(float)) {
			return rational_<T>(
				_CSTD cosf(static_cast<float>(x.numerator) / static_cast<float>(x.denominator))
				);
		} else {
			return rational_<T>(
				_CSTD cos(static_cast<double>(x.numerator) / static_cast<double>(x.denominator))
				);
		}
	}

	template<typename T> inline
	rational_<T> tan(rational_<T> x) {
		if constexpr (sizeof(T) <= sizeof(float)) {
			return rational_<T>(
				_CSTD tanf(static_cast<float>(x.numerator) / static_cast<float>(x.denominator))
				);
		} else {
			return rational_<T>(
				_CSTD tan(static_cast<double>(x.numerator) / static_cast<double>(x.denominator))
				);
		}
	}

	template<typename T> inline
	rational_<T> asin(rational_<T> x) {
		if constexpr (sizeof(T) <= sizeof(float)) {
			return rational_<T>(
				_CSTD asinf(static_cast<float>(x.numerator) / static_cast<float>(x.denominator))
				);
		} else {
			return rational_<T>(
				_CSTD asin(static_cast<double>(x.numerator) / static_cast<double>(x.denominator))
				);
		}
	}

	template<typename T> inline
	rational_<T> acos(rational_<T> x) {
		if constexpr (sizeof(T) <= sizeof(float)) {
			return rational_<T>(
				_CSTD acosf(static_cast<float>(x.numerator) / static_cast<float>(x.denominator))
				);
		} else {
			return rational_<T>(
				_CSTD acos(static_cast<double>(x.numerator) / static_cast<double>(x.denominator))
				);
		}
	}

	template<typename T> inline
	rational_<T> atan(rational_<T> x) {
		if constexpr (sizeof(T) <= sizeof(float)) {
			return rational_<T>(
				_CSTD atanf(static_cast<float>(x.numerator) / static_cast<float>(x.denominator))
				);
		} else {
			return rational_<T>(
				_CSTD atan(static_cast<double>(x.numerator) / static_cast<double>(x.denominator))
				);
		}
	}

}// namespace clmagic