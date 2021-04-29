#pragma once
/*{ "clmagic/calculation/fundamental/real":{
  "Author":"LongJiangnan",
  "Date":"2019-2021",
  "License":"Please identify Author",

  "
	decision-tree: {
	<!!! Only choose one, or make the same mistake>
		method1: { using template function,
			true: { 
				less code, 
				can still overrided 
			},
			false: {
				must replace all of it,
				consider very more usage enviroment,
				exist difficult and uncertain to optimize the special enviroment
					_Ty min(_Ty, _Ty)
					_Ty min(const _Ty&, const _Ty&)
					Which match?
			}
		},
		method2: { using override function,
			true: {
				more certain,
				can replace only part of it, not all of it
			},
			false: {
				more code,
				same meaning has to implement many times
			}
		}
  "
} }*/


#include <math.h>
#include "float.h"
#include "rational.h"
#include <numbers>
namespace calculation {
	using _CSTD abs;
	using _CSTD floor;
	inline float fract(float x) { return x - floor(x); }
	inline double fract(double x) { return x - floor(x); }
	using _CSTD ceil;
	using _CSTD round;

	using _CSTD pow;
	inline float pow(float x, int power) { return _CSTD powf(x, float(power)); }
	inline double pow(double x, int power) { return _CSTD pow(x, power); }
	inline int mod(int x, int y) { return x % y; }
	inline char mod(char x, char y) { return x % y; }
	inline float mod(float x, float y) { return fmodf(x, y); }
	inline double mod(double x, double y) { return fmod(x, y); }
	inline long long mod(long long x, long long y) { return x % y; }
	inline unsigned int mod(unsigned int x, unsigned int y) { return x % y; }
	inline unsigned char mod(unsigned char x, unsigned char y) { return x % y; }
	inline unsigned long long mod(unsigned long long x, unsigned long long y) { return x % y; }
	using _CSTD exp;
	using _CSTD log;
	using _CSTD sqrt;
	using _CSTD cbrt;
	
	using _CSTD sin;
	using _CSTD cos;
	using _CSTD tan;
	using _CSTD asin;
	using _CSTD acos;
	using _CSTD atan;

	/*<Theorem>
		_Ty a = any_value
		sin(a)*sin(a) + cos(a)*cos(a) = 1

		sin(a + 2Pi*k) = sin(a)
		cos(a + 2Pi*k) = cos(a)
		tan(a + 2Pi*k) = tan(a)
		cot(a + 2Pi*k) = cot(a)

		sin(a + Pi) = -sin(a)
		cos(a + Pi) = -cos(a)
		tan(a + Pi) =  tan(a)
		cot(a + Pi) =  cot(a)

		sin(-a) = -sin(a)
		cos(-a) =  cos(a)
		tan(-a) = -tan(a)
		cot(-a) = -cot(a)

		sin(a + Pi/2) =  cos(a)
		cos(a + Pi/2) = -sin(a)
		tan(a + Pi/2) = -cot(a)
		cot(a + Pi/2) = -tan(a)

		_Ty b = any_value
		sin(a + b) = sin(a)*cos(b) + cos(a)*sin(b)
		sin(a - b) = sin(a)*cos(b) - cos(a)*sin(b)
		cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)
		cos(a - b) = cos(a)*cos(b) + sin(a)*sin(b)
		tan(a + b) = ( tan(a)+tan(b) ) / ( 1 - tan(a)*tan(b) )
		tan(a - b) = ( tan(a)-tan(b) ) / ( 1 + tan(a)*tan(b) )
	</Theorem>*/

	/*<round-error>
		<example> cos(Pi/2) = -0.00000004F,
				  tan(angle / cos(Pi/2))
		</example>
		<avoid> tan2(angle, cos(Pi/2)), if(result < 0) result+=2Pi </avoid>
	<round-error>*/
	
#ifndef __CLMAGIC_CALULATION_SUM
#define __CLMAGIC_CALULATION_SUM
	template<typename Integer, typename Function>
	auto sum(Integer start, Integer end, Function f) -> decltype( f(start)+f(end) ) {
		auto result = f( start );
		for (Integer i = start+1; i <= end; ++i) {
			result += f( i );
		}

		return std::move( result );
	}
#endif

	template<typename Integer, typename Function>
	auto product(Integer start, Integer end, Function f) -> decltype( f(start)*f(end) ){
		auto result = f( start );
		for (Integer i = start+1; i <= end; ++i) {
			result *= f( i );
		}

		return std::move( result );
	}

	template<typename Real = double>
	Real factorial(uint32_t order) {
		Real result = 1;
		for (uint32_t i = 2; i <= order; ++i) {
			result *= i;
		}

		return result;
	}

	template<typename Real = double>
	Real approximately_factorial(uint32_t n) {// Mathematician:"Stirling"
		using std::numbers::pi;
		using std::numbers::e;
		return static_cast<Real>(
			sqrt(2*pi * n) * ::pow(n/e, n)
			);
	}


	/* minimum value */
	using std::min;
	
	/* maximum value */
	using std::max;
	
	/* value clamp to [lower, upper] */
	template<typename Ty1, typename Ty2> inline
	Ty1 clamp(const Ty1& value, const Ty2 lower, const Ty2 upper) {
		if (value < lower) {
			return static_cast<Ty1>(lower);
		}

		if (value > upper) {
			return static_cast<Ty1>(upper);
		}

		return value;
	}
	
	/* value from [lower,upper] to [new_lower, new_upper] */
	template<typename Ty1, typename Ty2, typename Ty3> inline
	Ty1 remap(const Ty1& value, const Ty2& lower, const Ty2& upper, const Ty3& new_lower, const Ty3& new_upper) {
		/* "no-boundary-error": { 
			"upper": "if value == lower, (lower - lower) / (upper - lower) * (new_upper - new_lower) + new_lower)
					                    = 0 + new_lower" 
			"lower": "if value == upper, (upper - lower) / (upper - lower) * (new_upper - new_lower) + new_lower)
										= 1 * (new_upper - new_lower) + new_lower"
		} */
		//assert(lower <= value && value <= upper);
		return static_cast<Ty1>((value - lower) / (upper - lower) * (new_upper - new_lower) + new_lower);
	}


	/* start to end by parameter(t) */
	template<typename Ty1, typename Ty2> inline
	Ty1 lerp(const Ty1& start, const Ty1& end, const Ty2& t) {
		return static_cast<Ty1>(start + (end - start) * t);
	}
	
	/* line(start0,start1) to line(end0,end1) by parameter(tX,tY) */
	template<typename Ty1, typename Ty2> inline
	Ty1 bilerp(const Ty1& start0, const Ty1& start1, const Ty1 end0, const Ty1 end1, const Ty2 tX, const Ty2 tY) {
		return static_cast<Ty1>(lerp( lerp(start0,start1,tX), lerp(end0, end1,tX), tY ));
	}
	
	/* plane(start00,start10,start01,start11) to plane(end00,end10,end01,end11) by parameter(tX,tY,tZ) */
	template<typename Ty1, typename Ty2> inline
	Ty1 trilerp(const Ty1& start00, const Ty1& start10, const Ty1& start01, const Ty1 start11, const Ty1& end00, const Ty1& end10, const Ty1& end01, const Ty1 end11, 
		const Ty2 tX, const Ty2 tY, const Ty2 tZ) {
		return static_cast<Ty1>(lerp( bilerp(start00,start10,start01,start11,tX,tY), bilerp(end00,end10,end01,end11,tX,tY), tZ ));
	}
	
	// bilinear sample function(f) at (u,v)
	template<typename Fn, typename Ty> requires requires(Fn __f, Ty __t) { __f(__t, __t); }
	auto bilersmp(const Fn& latticef, Ty u, Ty v) {
		// setup lerp parameters
		size_t x = static_cast<size_t>(floor(u));
		size_t y = static_cast<size_t>(floor(v));
		Ty tX = fract(u);
		Ty tY = fract(v);
			
		// lerp four lattices
		auto sample = latticef(x,y);
		if (tX != 0) {
			sample = lerp(sample, latticef(x+1,y), tX);
			if (tY != 0) {
				sample = lerp(sample, lerp(latticef(x,y+1),latticef(x+1,y+1),tX), tY);
			}
		} else {
			if (tY != 0) {
				sample = lerp(sample, latticef(x,y+1), tY);
			}
		}

		return sample;
	}
	
	template<typename Fn, typename Ty> requires requires(Fn __f, Ty __t) { __f(__t, __t, __t); }
	auto trilersmp(const Fn& latticef, Ty u, Ty v, Ty w) {
		// setup lerp parameters
		size_t x = static_cast<size_t>(floor(u));
		size_t y = static_cast<size_t>(floor(v));
		size_t z = static_cast<size_t>(floor(w));
		Ty tX = fract(u);
		Ty tY = fract(v);
		Ty tZ = fract(w);
			
		// lerp eight lattices
		auto sample = latticef(x,y,z);
		if (tX != 0) {
			size_t x1 = x + 1;
			sample = lerp(sample, latticef(x1,y,z), tX);
			if (tY != 0) {
				size_t y1 = y + 1;
				sample = lerp(sample, lerp(latticef(x,y1,z),latticef(x1,y1,z),tX), tY);
				if (tZ != 0) {
					size_t z1 = z + 1;
					sample = lerp(sample, bilerp(latticef(x,y,z1),latticef(x1,y,z1),latticef(x,y1,z1),latticef(x1,y1,z1),tX,tY), tZ);
				}
			} else {
				if (tZ != 0) {
					// lerp XZ
					size_t z1 = z + 1;
					sample = lerp(sample, lerp(latticef(x,y,z1),latticef(x1,y,z1),tX), tZ);
				}
			}
		} else {

			if (tY != 0) {
				size_t y1 = y + 1;
				sample = lerp(sample, latticef(x,y1,z), tY);
				if (tZ != 0) {
					size_t z1 = z + 1;
					sample = lerp(sample, lerp(latticef(x,y,z1),latticef(x,y1,z1),tY), tZ);
				}
			} else {
				if (tZ != 0) {
					sample = lerp(sample, latticef(x,y,z+1), tZ);
				}
			}

		}

		return sample;
	}



	/* wrapping a lattice-function
	* latticef(domain:[any]) -> range:[any]
	* normalized_latticef(domain:[0,1]x[0,1]) -> range:[any]
	*/
	template<typename Lattice> requires requires(Lattice __f) { __f(size_t(), size_t()); }
	class NormalizedSampler2D {
	public:
		using lattice_type = Lattice;
		using result_type = decltype(Lattice()(size_t(),size_t()));
		
		const lattice_type& lattice;
		size_t lattice_rows;
		size_t lattice_cols;

		template<typename Ty>
		result_type operator()(Ty u, Ty v) const {
			assert( 0 <= u && u <= 1 );
			assert( 0 <= v && v <= 1 );
			return bilersmp(lattice, u*(lattice_cols-1), v*(lattice_rows-1));
		}

	public:
		NormalizedSampler2D(const lattice_type& _lattice, size_t _rows, size_t _cols)
			: lattice(_lattice), lattice_rows(_rows), lattice_cols(_cols) {}
		size_t rows() const { return lattice_rows; }
		size_t cols() const { return lattice_cols; }
	};

	template<typename Lattice> requires requires(Lattice __f) { __f(size_t(), size_t(), size_t()); }
	class NormalizedSampler3D {
	public:
		using lattice_type = Lattice;
		using result_type = decltype(Lattice()(size_t(),size_t(), size_t()));
		
		const lattice_type& lattice;
		size_t lattice_rows;
		size_t lattice_cols;
		size_t lattice_slices;

		template<typename Ty>
		result_type operator()(Ty u, Ty v, Ty w) const {
			assert( 0 <= u && u <= 1 );
			assert( 0 <= v && v <= 1 );
			assert( 0 <= w && w <= 1 );
			return trilersmp(lattice, u*(lattice_cols-1), v*(lattice_rows-1), w*(lattice_slices-1));
		}

	public:
		NormalizedSampler3D(const lattice_type& _lattice, size_t _rows, size_t _cols, size_t _slices)
			: lattice(_lattice), lattice_rows(_rows), lattice_cols(_cols), lattice_slices(_slices) {}
		size_t rows() const { return lattice_rows; }
		size_t cols() const { return lattice_cols; }
		size_t slices() const { return lattice_slices; }
	};

	/* wrapping a lattice-function
	* latticef(domain:[any]) -> range:[any]
	* heightmap(domain:[horizontal_lowest,horizontal_max]) -> range:[vertical_lowest,vertical_max]
	*/
	template<typename Lattice, typename Length = double> requires requires(Lattice __f) { __f(size_t(), size_t()); }
	class HeightmapSampler {
	public:
		using lattice_type = Lattice;
		using value_type = decltype(Lattice()(size_t(),size_t()));
		using length_type = Length;

		const lattice_type& lattice;
		size_t lattice_row_backindex;
		size_t lattice_col_backindex;
		value_type lattice_value_lowest;
		value_type lattice_value_max;

		length_type horizontal_lowest[2];
		length_type horizontal_max[2];
		length_type vertical_lowest;
		length_type vertical_max;

		template<typename Ty>
		length_type operator()(const Ty& x, const Ty& z) const {
			assert( horizontal_lowest[0] <= x && x <= horizontal_max[0] );
			assert( horizontal_lowest[1] <= z && z <= horizontal_max[1] );

			Ty u = remap(x, horizontal_lowest[0],horizontal_max[0],
				size_t(0),lattice_row_backindex);

			Ty v = remap(z, horizontal_lowest[1],horizontal_max[1],
				size_t(0),lattice_col_backindex);

			return remap(bilersmp(lattice, u, v), lattice_value_lowest,lattice_value_max,
				vertical_lowest,vertical_max);
		}

	public:
		HeightmapSampler(const lattice_type& _lattice, size_t _rows, size_t _cols, value_type _lowest, value_type _max,
			length_type arg0_lowest, length_type arg0_max,
			length_type arg1_lowest, length_type arg1_max,
			length_type result_lowest, length_type result_max) : lattice(_lattice) {
			lattice_row_backindex = _rows - 1;
			lattice_col_backindex = _cols - 1;
			lattice_value_lowest = _lowest;
			lattice_value_max = _max;
			horizontal_lowest[0] = arg0_lowest;
			horizontal_max[0] = arg0_max;
			horizontal_lowest[1] = arg1_lowest;
			horizontal_max[1] = arg1_max;
			vertical_lowest = result_lowest;
			vertical_max = result_max;
		}

		HeightmapSampler(const lattice_type& _lattice, size_t _rows, size_t _cols,
			length_type arg0_lowest, length_type arg0_max,
			length_type arg1_lowest, length_type arg1_max,
			length_type result_lowest, length_type result_max)
			: HeightmapSampler(_lattice, _rows, _cols, 0, 1, arg0_lowest, arg0_max, arg1_lowest, arg1_max, result_lowest, result_max) {}
		
		size_t rows() const { 
			return lattice_row_backindex + 1;
		}
		
		size_t cols() const { 
			return lattice_col_backindex + 1;
		}
	};
	

	template<typename Ty> 
	struct Range {
		Ty lower;
		Ty upper;
		Range() = default;
		Range(Ty _lower, Ty _upper) : lower(_lower), upper(_upper) {}

		template<typename Ty2>
		bool operator()(const Ty2& value) const {
			return lower <= value && value <= upper;
		}
	};
	
	template<typename Ty1, typename Ty2> inline
	Ty1 clamp(const Ty1& value, const Range<Ty2>& range) {
		return clamp(value, range.lower, range.upper);
	}
	
	template<typename Ty1, typename Ty2, typename Ty3> inline
	Ty1 remap(const Ty1& value, const Range<Ty2>& input_range, const Range<Ty3>& output_range) {
		assert( input_range(value) == true );
		return remap(value, input_range.lower, input_range.upper, output_range.lower, output_range.upper);
	}
}