#pragma once
/*{ "clmagic/calculation/fundamental/float":{
  "Author": "LongJiangnan",
  "Date": "2019-2021",
  "License": "Please identify Author",
} }*/


#include <initializer_list>
#include <iterator>
#include <functional>
namespace calculation {
	template<typename InputIterator, typename Array>
	void set_array(InputIterator psource, Array& destination) {
		auto pdest = destination.begin();
		auto pdest_end = destination.end();
		for ( ; pdest != pdest_end; ++pdest, ++psource) {
			*pdest = *psource;
		}
	}

	// transform X[i] into Y[i] using f
	template<typename InArray, typename OutArray, typename UnaryOp>
	OutArray& transform_static_array(const InArray& X, OutArray& Y, UnaryOp f) {
		assert(X.size() >= OutArray::size());

		if constexpr (OutArray::size() == 1) {
			Y[0] = f(X[0]);
		} else if constexpr (OutArray::size() == 2) {
			Y[0] = f(X[0]); Y[1] = f(X[1]);
		} else if constexpr (OutArray::size() == 3) {
			Y[0] = f(X[0]); Y[1] = f(X[1]); Y[2] = f(X[2]);
		} else if constexpr (OutArray::size() == 4) {
			Y[0] = f(X[0]); Y[1] = f(X[1]); Y[2] = f(X[2]); Y[3] = f(X[3]);
		} else {
			auto first = X.begin();
			auto last = X.end();
			auto dest = Y.begin();
			for ( ; first != last; ++first, ++dest) {
				*dest = f(*first);
			}
		}

		return Y;
	}
	// transform X[i] and Y[i] into Z[i] using f
	template<typename InArray1, typename InArray2, typename OutArray, typename BinaryOp>
	OutArray& transform_static_array(const InArray1& X, const InArray2& Y, OutArray& Z, BinaryOp f) {
		assert(X.size() >= OutArray::size());
		assert(Y.size() >= OutArray::size());

		if constexpr (OutArray::size() == 1) {
			Z[0] = f(X[0],Y[0]);
		} else if constexpr (OutArray::size() == 2) {
			Z[0] = f(X[0],Y[0]); Z[1] = f(X[1],Y[1]);
		} else if constexpr (OutArray::size() == 3) {
			Z[0] = f(X[0],Y[0]); Z[1] = f(X[1],Y[1]); Z[2] = f(X[2],Y[2]);
		} else if constexpr (OutArray::size() == 4) {
			Z[0] = f(X[0],Y[0]); Z[1] = f(X[1],Y[1]); Z[2] = f(X[2],Y[2]); Z[3] = f(X[3],Y[3]);
		} else {
			auto first1 = X.begin();
			auto last1  = X.end();
			auto first2 = Y.begin();
			auto dest   = Z.begin();
			for ( ; first1 != last1; ++first1, ++first2, ++dest) {
				*dest = f(*first1, *first2);
			}
		}

		return Z;
	}

	// f(X[i])
	template<typename OutArray, typename InArray, typename UnaryOp>
	OutArray make_static_array(const InArray& X, UnaryOp f) {
		assert(X.size() >= OutArray::size());

		if constexpr (OutArray::size() == 1) {
			return OutArray{ f(X[0]) };
		} else if constexpr (OutArray::size() == 2) {
			return OutArray{ f(X[0]), f(X[1]) };
		} else if constexpr (OutArray::size() == 3) {
			return OutArray{ f(X[0]), f(X[1]), f(X[2]) };
		} else if constexpr (OutArray::size() == 4) {
			return OutArray{ f(X[0]), f(X[1]),  f(X[2]),  f(X[3]) };
		} else {
			OutArray Y;
			auto first = X.begin();
			auto last = X.end();
			auto dest = Y.begin();
			for ( ; first != last; ++first, ++dest) {
				*dest = f(*first);
			}
			return Y;
		}
	}
	// f(X[i],Y[i])
	template<typename OutArray, typename InArray1, typename InArray2, typename BinaryOp>
	OutArray make_static_array(const InArray1& X, const InArray2& Y, BinaryOp f) {
		assert(X.size() >= OutArray::size());
		assert(Y.size() >= OutArray::size());

		if constexpr (OutArray::size() == 1) {
			return OutArray{ f(X[0],Y[0]) };
		} else if constexpr (OutArray::size() == 2) {
			return OutArray{ f(X[0],Y[0]), f(X[1],Y[1]) };
		} else if constexpr (OutArray::size() == 3) {
			return OutArray{ f(X[0],Y[0]), f(X[1],Y[1]), f(X[2],Y[2]) };
		} else if constexpr (OutArray::size() == 4) {
			return OutArray{ f(X[0], Y[0]), f(X[1],Y[1]), f(X[2],Y[2]), f(X[3],Y[3]) };
		} else {
			OutArray Z;
			auto first1 = X.begin();
			auto last1  = X.end();
			auto first2 = Y.begin();
			auto dest   = Z.begin();
			for ( ; first1 != last1; ++first1, ++first2, ++dest) {
				*dest = f(*first1, *first2);
			}
			return Z;
		}
	}


	// static array
	template<typename T, size_t N, typename ExPo = void>
	struct arrayX {
	public:
		using array_type = arrayX<T,N,ExPo>;
		using scalar_type = T;
		using execution_policy = ExPo;

		static constexpr 
		size_t size() {
			return N;
		}

		scalar_type _Mydata[size()];

		const scalar_type* data() const { 
			return _Mydata; 
		}
		scalar_type* data() { 
			return _Mydata; 
		}
		
		const scalar_type* begin() const { 
			return _Mydata; 
		}
		scalar_type* begin() { 
			return _Mydata; 
		}
		
		const scalar_type* end() const { 
			return _Mydata + size(); 
		}
		scalar_type* end() { 
			return _Mydata + size(); 
		}
		
		const scalar_type& at(size_t i) const { 
			return _Mydata[i]; 
		}
		scalar_type& at(size_t i) { 
			return _Mydata[i]; 
		}
		
		const scalar_type& operator[](size_t i) const { 
			return _Mydata[i]; 
		}
		scalar_type& operator[](size_t i) { 
			return _Mydata[i]; 
		}

	public:
		template<typename rT, size_t rN, typename rExPo> explicit 
		operator arrayX<rT,rN,rExPo>() const {
			arrayX<rT,rN,rExPo> result;
			std::copy(_Mydata, _Mydata + std::min(size(),result.size()), 
				result._Mydata);
			return std::move( result );
		}

		array_type operator-() const {
			// make static array by [*this] with std::negate<void>()
			return make_static_array<array_type>(*this, std::negate<>());
		}
		
		array_type& operator+=(const array_type& other) {
			transform_static_array(*this, other, *this, std::plus<>());
			return *this;
		}
		array_type operator+(const array_type& other) const {
			// make static array from [*this] and [other] with std::plus<void>()
			return make_static_array<array_type>(*this, other, std::plus<>());
		}
		
		array_type& operator-=(const array_type& other) {
			transform_static_array(*this, other, *this, std::minus<>());
			return *this;
		}
		array_type operator-(const array_type& other) const {
			// make static array from [*this] and [other] with std::minus<void>()
			return make_static_array<array_type>(*this, other, std::minus<>());
		}
		
		array_type& operator*=(const array_type& other) {
			transform_static_array(*this, other, (*this), std::multiplies<>());
			return *this;
		}
		array_type operator*(const array_type& other) const {
			// make static array from [*this] and [other] with std::multiplies<void>()
			return make_static_array<array_type>(*this, other, std::multiplies<>());
		}
		
		array_type& operator/=(const array_type& other) {
			transform_static_array(*this, other, (*this), std::divides<>());
			return *this;
		}
		array_type operator/(const array_type& other) const {
			// make static array from [*this] and [other] with std::divides<void>()
			return make_static_array<array_type>(*this, other, std::divides<>());
		}
		
		array_type& operator%=(const array_type& other) {
			if constexpr (std::is_same_v<scalar_type, float>) {
				transform_static_array(*this, other, (*this), [](float x, float y){ return _CSTD fmodf(x, y); });
			} else if constexpr (std::is_same_v<scalar_type, double>) {
				transform_static_array(*this, other, (*this), [](double x, double y){ return _CSTD fmod(x, y); });
			} else {
				transform_static_array(*this, other, (*this), std::modulus<>());
			}
			return *this;
		}
		array_type operator%(const array_type& other) const {
			// make static array from [*this] and [other] with std::modulus<void>()
			if constexpr (std::is_same_v<scalar_type, float>) {
				return make_static_array<array_type>(*this, other, [](float x, float y) { return _CSTD fmodf(x, y); });
			} else if constexpr (std::is_same_v<scalar_type, double>) {
				return make_static_array<array_type>(*this, other, [](double x, double y) { return _CSTD fmod(x, y); });
			} else {
				return make_static_array<array_type>(*this, other, std::modulus<>());
			}
		}

		array_type& operator+=(const scalar_type& adder) {
			transform_static_array( *this, *this, [&adder](const scalar_type& x){ 
				return x + adder; } );
			return *this;
		}
		array_type operator+(const scalar_type& adder) const {
			return make_static_array<array_type>( *this, [&adder](const scalar_type& x){ 
				return x + adder; } );
		}

		array_type& operator-=(const scalar_type& substractor) {
			transform_static_array( *this, *this, [&substractor](const scalar_type& x){
				return x - substractor; } );
			return *this;
		}
		array_type operator-(const scalar_type& substractor) const {
			return make_static_array<array_type>( *this, [&substractor](const scalar_type& x){
				return x - substractor; } );
		}

		array_type& operator*=(const scalar_type& multiplier) {
			transform_static_array( *this, *this, [&multiplier](const scalar_type& x){
				return x * multiplier; } );
			return *this;
		}
		array_type operator*(const scalar_type& multiplier) const {
			return make_static_array<array_type>( *this, [&multiplier](const scalar_type& x){
				return x * multiplier; } );
		}

		array_type& operator/=(const scalar_type& divisor) {
			transform_static_array( *this, *this, [&divisor](const scalar_type& x){
				return x / divisor; } );
			return *this;
		}
		array_type operator/(const scalar_type& divisor) const {
			return make_static_array<array_type>( *this, [&divisor](const scalar_type& x){
				return x / divisor; } );
		}

		array_type& operator%=(const scalar_type& divisor) {
			if constexpr ( std::is_same_v<scalar_type, float> ) {
				transform_static_array( *this, *this, [&divisor](const scalar_type& x){
					return _CSTD fmodf(x, divisor); } );
			} else if constexpr ( std::is_same_v<scalar_type, double> ) {
				transform_static_array( *this, *this, [&divisor](const scalar_type& x){
					return _CSTD fmod(x, divisor); } );
			} else {			
				transform_static_array( *this, *this, [&divisor](const scalar_type& x){
					return x % divisor; } );
			}
			return *this;
		}
		array_type operator%(const scalar_type& divisor) const {
			if constexpr ( std::is_same_v<scalar_type, float> ) {
				return make_static_array<array_type>( *this, [&divisor](const scalar_type& x){
					return _CSTD fmodf(x, divisor); } );
			} else if constexpr ( std::is_same_v<scalar_type, double> ) {
				return make_static_array<array_type>( *this, [&divisor](const scalar_type& x){
					return _CSTD fmod(x, divisor); } );
			} else {
				return make_static_array<array_type>( *this, [&divisor](const scalar_type& x){
					return x % divisor; } );
			}
		}

		array_type operator==(const array_type& other) const {
			return make_static_array<array_type>( *this, other, [](const scalar_type& x, const scalar_type& y) {
				return static_cast<scalar_type>(x == y); 
				} );
		}
		array_type operator==(const scalar_type& target) const {
			return make_static_array<array_type>( *this, [&target](const scalar_type& x) {
				return static_cast<scalar_type>(x == target);
				} );
		}

		array_type operator!=(const array_type& other) const {
			return make_static_array<array_type>( *this, other, [](const scalar_type& x, const scalar_type& y) {
				return static_cast<scalar_type>(x != y);
				} );
		}
		array_type operator!=(const scalar_type& target) const {
			return make_static_array<array_type>( *this, [&target](const scalar_type& x) {
				return static_cast<scalar_type>(x != target);
				} );
		}

		array_type operator<(const array_type& other) const {
			return make_static_array<array_type>( *this, other, [](const scalar_type& x, const scalar_type& y) {
				return static_cast<scalar_type>(x < y);
				} );
		}
		array_type operator<(const scalar_type& upper) const {
			return make_static_array<array_type>( *this, [&upper](const scalar_type& x) {
				return static_cast<scalar_type>(x < upper);
				} );
		}

		array_type operator>(const array_type& other) const {
			return make_static_array<array_type>( *this, other, [](const scalar_type& x, const scalar_type& y) {
				return static_cast<scalar_type>(x > y);
				} );
		}
		array_type operator>(const scalar_type& lower) const {
			return make_static_array<array_type>( *this, [&lower](const scalar_type& x) {
				return static_cast<scalar_type>(x > lower);
				} );
		}

		array_type operator<=(const array_type& other) const {
			return make_static_array<array_type>( *this, other, [](const scalar_type& x, const scalar_type& y) {
				return static_cast<scalar_type>(x <= y);
				} );
		}
		array_type operator<=(const scalar_type& upper) const {
			return make_static_array<array_type>( *this, [&upper](const scalar_type& x) {
				return static_cast<scalar_type>(x <= upper);
				} );
		}

		array_type operator>=(const array_type& other) const {
			return make_static_array<array_type>( *this, other, [](const scalar_type& x, const scalar_type& y) {
				return static_cast<scalar_type>(x >= y);
				} );
		}
		array_type operator>=(const scalar_type& lower) const {
			return make_static_array<array_type>( *this, [&lower](const scalar_type& x) {
				return static_cast<scalar_type>(x >= lower);
				} );
		}

		friend array_type operator+(const scalar_type& adder, const array_type& the_array) {
			return the_array + adder;
		}
		friend array_type operator*(const scalar_type& adder, const array_type& the_array) {
			return the_array * adder;
		}
	};

	template<typename T, size_t N> inline
	arrayX<T, N> min(const arrayX<T,N>& X, const arrayX<T,N>& Y) {
		return make_static_array< arrayX<T,N> >( X, Y, [](const T& x, const T& y){ 
			return min(x,y); } );
	}

	template<typename T, size_t N> inline
	arrayX<T, N> max(const arrayX<T,N>& X, const arrayX<T,N>& Y) {
		return make_static_array< arrayX<T,N> >( X, Y, [](const T& x, const T& y){ 
			return max(x,y); } );
	}

	template<typename T, size_t N> inline
	arrayX<T,N> abs(const arrayX<T,N>& X) {
		return make_static_array<arrayX<T,N>>( X, [](const T& x){ 
			return abs(x); } );
	}
	
	template<typename T, size_t N> inline
	arrayX<T,N> floor(const arrayX<T,N>& X) {
		return make_static_array<arrayX<T,N>>( X, [](const T& a){ 
			return floor(a); } );
	}

	template<typename T, size_t N> inline
	arrayX<T,N> fract(const arrayX<T,N>& X) {
		return make_static_array<arrayX<T,N>>( X, [](const T& a){ 
			return fract(a); } );
	}

	template<typename T, size_t N> inline
	arrayX<T,N> ceil(const arrayX<T,N>& X) {
		return make_static_array<arrayX<T,N>>( X, [](const T& a){ 
			return ceil(a); } );
	}

	template<typename T, size_t N> inline
	arrayX<T,N> round(const arrayX<T,N>& X) {
		return make_static_array<arrayX<T,N>>( X, [](const T& a){ 
			return round(a); } );
	}
	
	template<typename T, size_t N> inline
	arrayX<T,N> sqrt(const arrayX<T,N>& X) {
		return make_static_array<arrayX<T, N>>( X, [](const T& x){ 
			return sqrt(x); } );
	}
	
	template<typename T, size_t N> inline
	arrayX<T,N> cbrt(const arrayX<T,N>& X) {
		return make_static_array<arrayX<T, N>>( X, [](const T& x){ 
			return cbrt(x); } );
	}

	template<typename T, size_t N> inline
	arrayX<T,N> exp(const arrayX<T,N>& X) {
		return make_static_array<arrayX<T, N>>( X, [](const T& x){ 
			return exp(x); } );
	}

	template<typename T, size_t N> inline
	auto dot(const arrayX<T,N>& v1, const arrayX<T,N>& v2) -> typename arrayX<T,N>::scalar_type {
		if constexpr (N == 2) {
			return v1[0]*v2[0] + v1[1]*v2[1];
		} else if constexpr (N == 3) {
			return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
		} else if constexpr (N == 4) {
			return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3];
		} else {
			typename arrayX<T,N>::scalar_type sigma = v1[0] * v2[0];
			for (size_t i = 1; i != N; ++i) {
				sigma += v1[i] * v2[i];
			}

			return sigma;
		}
	}


	template<typename T>
		using array2 = arrayX<T, 2>;
	template<typename T>
		using array3 = arrayX<T, 3>;
	template<typename T>
		using array4 = arrayX<T, 4>;
	template<typename T>
		using array5 = arrayX<T, 5>;
	template<typename T>
		using array6 = arrayX<T, 6>;
	template<typename T>
		using array7 = arrayX<T, 7>;
	template<typename T>
		using array8 = arrayX<T, 8>;
	using array2d = arrayX<double, 2>;
	using array3d = arrayX<double, 3>;
	using array4d = arrayX<double, 4>;
	using array5d = arrayX<double, 5>;
	using array6d = arrayX<double, 6>;
	using array7d = arrayX<double, 7>;
	using array8d = arrayX<double, 8>;
}// namespace array