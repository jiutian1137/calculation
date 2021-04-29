#pragma once
/*{ "clmagic/calculation/formulas/legendre":{

"Mathematician":[
  "Adrien Marie Legendre",
  "Francesco Giacomo Tricomi",
  "Frank G. Lether",
  "Gatteschi, L."
] ,

"Reference":[
{ "Url":"https://www.physics.uoguelph.ca/chapter-3-legendre-polynomials",
  "Author":"?",
  "Date":"?",
  "Desc":"understand Legendre polynomials" } ,

{ "Paper":"On the construction of Gauss-Legendre quadrature rules" ,
  "Author":"Frank G. Lether",
  "Date":"1978",
  "Desc":"5th iteration method" } ,

{ "Lib":"boost::math::quadrature::gauss<T>",
  "Author":"John Maddock",
  "Date":"2015" } ,

{ "Lib":"boost::math::legendre_*",
  "Author":"John Maddock",
  "Date":"2006" }
] ,

"Author":"LongJiangnan" ,

"License":"Please Identify Mathematician"

} }*/


#include <cassert>
#include <numeric>
#if __has_include("clmagic/calculation/fundamental/real.h") == 1
#include "../fundamental/real.h"
#else
#include <cmath>
#include <utility>
namespace calculation {
	using std::min;
	using std::max;
}
#endif
namespace calculation {
	using LegendrePolynomialOrder = int;
	
	/*Legendre Polynomial First Kind
					   1             d(L)th( pow(x*x-1, n) ) 
	* P(L) = ------------------- * -------------------------
		   	   pow(2,L) * fact(L)         dx(L)th
	* P(0) = 1
	* P(1) = x
	* P(2) = (1/2) * (3*x^2 - 1)
	* P(3) = (1/2) * (5*x^3 - 3*x)
	* P(4) = (1/8) * (35*x^4 - 30*x^2 + 3)
	* P(5) = (1/8) * (63*x^5 - 70*x^3 + 15*x)
	* P(6) = (1/16) * (231*x^6 - 315*x^4 + 105*x - 5)
	   ...
	* Function Object
	*/ 
	struct LegendrePolynomial {
		LegendrePolynomialOrder order;
		
		template<typename Real>
		Real operator()(Real x)const {
			// evaluate P(x), recursion: Pnext = ( Pcurrent*x*(2*i+1) - Pprev*i ) / (i+1) }
			if (order == 0) {
				return static_cast<Real>(1);
			}

			Real Pprev = static_cast<Real>(1);
			Real Pcurrent = x;
			for (LegendrePolynomialOrder i = 1; i < order; ++i) {
				Real Pnext = ( Pcurrent * ( x * (2*i+1) ) - Pprev * i ) / (i+1);
				Pprev = Pcurrent;
				Pcurrent = Pnext;
			}
			return Pcurrent;
		}

		LegendrePolynomial() = default;
		LegendrePolynomial(LegendrePolynomialOrder L) : order(L) {}
		LegendrePolynomial& operator+=(int diff) {
			order += diff;
			return *this;
		}
		LegendrePolynomial& operator-=(int diff) {
			order -= diff;
			return *this;
		}
		LegendrePolynomial operator+(int diff) const {
			return LegendrePolynomial(*this) += diff;
		}
		LegendrePolynomial operator-(int diff) const {
			return LegendrePolynomial(*this) -= diff;
		}
		LegendrePolynomial operator-() const {
			return LegendrePolynomial(-order);
		}
	};

	/*d(P(x))/dx
	* differentiate P(x), recursion: dF(Pnext) = Pcurrent*(i+1) + dF(Pcurrent)*x
	*/
	template<typename Real>
	Real differentiate(LegendrePolynomial P, Real x, Real* evaluated = nullptr) {
		if (P.order == 0) {
			if (evaluated) { 
				*evaluated = static_cast<Real>(1); 
			}
			return static_cast<Real>(0);
		}

		Real Pprev = static_cast<Real>(1);
		Real Pcurrent = x;
		Real dPcurrent = static_cast<Real>(1);// derivative of x
		for (LegendrePolynomialOrder i = 1; i < P.order; ++i) {
			Real Pnext = ( Pcurrent * ( x * (2*i+1) ) - Pprev * i ) / (i+1);
			Real dPnext = Pcurrent * (i + 1) + dPcurrent * x;
			Pprev = Pcurrent;
			Pcurrent = Pnext;
			dPcurrent = dPnext;
		}

		if (evaluated) {
			*evaluated = Pcurrent;
		}
		return dPcurrent;
	}
	
	/*solve(P(x) = 0)
	* root_array_size is size_t(P.order)
	*/
	template<typename Real, std::enable_if_t<sizeof(Real) >= 8, int> = 0>
	size_t solve_zeros(LegendrePolynomial P, Real* root_array) {
		if ( P.order < 0 ) {
			return solve_zeros(LegendrePolynomial(-P.order - 1), root_array);
		}
		
		assert( root_array != nullptr );		
		const Real pi = static_cast<Real>(3.141592653589793);
		const Real eps = std::numeric_limits<Real>::epsilon();
		const uint32_t n = static_cast<uint32_t>(P.order);
		assert( n <= sqrt(std::numeric_limits<uint32_t>::max() / sqrt(384.0)) );// n <= 14491, uint32_t over flow
		uint32_t k = 1;

		for ( ; k <= 1; ++k) {
			// root = x + O(n^(-7)), Mathematician:"Gatteschi, L."
			Real jk = static_cast<Real>(2.404825557695772768622);
			Real a  = static_cast<Real>(1)/3 + n*n + n;
			Real x  = cos( jk/sqrt(a) * ( 1 - (jk*jk - 2)/(a*a * 360) ) );
			// 5-th order iterative method, Mathematician:"Frank G. Lether"
			Real xprev;
			do {
				xprev = x;
				Real Pnm1 = (P - 1)( x );
				Real Pn = P( x );
				Real v = Pn / (n*Pnm1 - n*x*Pn);
				Real B = static_cast<Real>(3 + n*(n+1))*(x*x)/3 - static_cast<Real>(1 + n*(n+1))/3;
				Real C = static_cast<Real>(6 + 5*n*(n+1))*(x*x*x)/6 - static_cast<Real>(4 + 5*n*(n+1))/6;
				x = x - (1-x)*(1+x)*v*(1 + v*(x + v*(B + C*v)));
			} while ( abs(x-xprev) >= eps*max(abs(x),eps) );
			// output
			root_array[k-1] = x;
		}

		for ( ; k <= n/3; ++k) {
			// root = x + O(n^(-5)), Mathematician:"Tricomi"
			Real v  = pi*(4*k-1) / (4*n+2);
			Real sv = sin(v);
			Real x  = cos(v) * ( (1 - static_cast<Real>(1)/(8U*n*n) + static_cast<Real>(1)/(8ULL*n*n*n)) - (39 - 28/(sv*sv))/(384ULL*n*n*n*n) );
			// 5-th order iterative method, Mathematician:"Frank G. Lether"
			Real xprev;
			do {
				xprev = x;
				Real Pnm1 = (P - 1)( x );
				Real Pn = P( x );
				Real v = Pn / (n*Pnm1 - n*x*Pn);
				Real B = static_cast<Real>(3 + n*(n+1))*(x*x)/3 - static_cast<Real>(1 + n*(n+1))/3;
				Real C = static_cast<Real>(6 + 5*n*(n+1))*(x*x*x)/6 - static_cast<Real>(4 + 5*n*(n+1))/6;
				x = x - (1-x)*(1+x)*v*(1 + v*(x + v*(B + C*v)));
			} while ( abs(x-xprev) >= eps*max(abs(x),eps) );
			// output
			root_array[k-1] = x;
		}

		for ( ; k <= n/2; ++k) {
			// root = x + O(n^(-4)), Mathematician:"Tricomi"
			Real v = pi*(4*k-1) / (4*n+2);
			Real x = cos(v) * (1 - static_cast<Real>(1)/(8U*n*n) + static_cast<Real>(1)/(8ULL*n*n*n));
			// 5-th order iterative method, Mathematician:"Frank G. Lether"
			Real xprev;
			do {
				xprev = x;
				Real Pnm1 = (P - 1)( x );
				Real Pn = P( x );
				Real v = Pn / (n*Pnm1 - n*x*Pn);
				Real B = static_cast<Real>(3 + n*(n+1))*(x*x)/3 - static_cast<Real>(1 + n*(n+1))/3;
				Real C = static_cast<Real>(6 + 5*n*(n+1))*(x*x*x)/6 - static_cast<Real>(4 + 5*n*(n+1))/6;
				x = x - (1-x)*(1+x)*v*(1 + v*(x + v*(B + C*v)));
			} while ( abs(x-xprev) >= eps*max(abs(x),eps) );
			// output
			root_array[k-1] = x;
		}

		if (n & 1) {// symmetry middle
			assert(k == n/2+1);
			root_array[k-1] = static_cast<Real>(0);
			++k;
		}

		for ( ; k <= n; ++k) {// symmetry
			root_array[k - 1] = -root_array[n+1-k - 1];
		}

		return static_cast<size_t>(P.order);
	}

	/*gauss_quadrature weight_array
	*/
	template<typename Real>
	void legendre_gauss_integrate_weights(LegendrePolynomial P, const Real* root_array, Real* weight_array) {
		assert( root_array != nullptr );
		assert( weight_array != nullptr );
		for (int i = 0, n = static_cast<int>(P.order); i < n; ++i) {
			Real x  = root_array[i];
			Real dP = differentiate(P, x);
			weight_array[i] = 2 / ( (1 - x*x) * dP*dP );
		}
	}

	/*gauss_quadrature
	* sum( (scale * wi) * f(center + ri * scale) ), Legendre polynomial zeros are symmetry
	*/
	template<typename F, typename Real>
	Real legendre_gauss_integrate(F f, Real a, Real b, const Real* root_array, const Real* weight_array, size_t array_size) {
		Real scale = (b - a) / 2;
		Real center = (a + b) / 2;
		size_t nhalf_i = array_size / 2;

		Real result = Real(0);
		if (array_size & 1) {
			Real xi = root_array[nhalf_i];
			Real wi = weight_array[nhalf_i];
			result += wi * f(center + xi * scale)/* *scale */; ++nhalf_i;
		}
		for ( ; nhalf_i != array_size; ++nhalf_i) {
			Real xi = root_array[nhalf_i];
			Real wi = weight_array[nhalf_i];
			result += wi * ( f(center + xi * scale)/* *scale */ +
							 f(center - xi * scale)/* *scale */ );
		}

		return result * scale;
	}



	struct LegendrePolynomial2 {
		LegendrePolynomialOrder order;
		LegendrePolynomial2(LegendrePolynomialOrder L) : order(L) {}
	};

	// ...
}

/* GaussQuadratureParameter:{
LegendrePolynomial(2):
	double x0 = -0.57735026918963; // -sqrt(3.0)/3
	double x1 =  0.57735026918963; //  sqrt(3.0)/3
	double w0 =  1.00000000000000;
	double w1 =  1.00000000000000;

LegendrePolynomial(3):
	double x0 = -0.77459666924148; // -sqrt(15.0)/5
	double x1 =  0.00000000000000; //  0
	double x2 =  0.77459666924148; //  sqrt(15.0)/5
	double w0 =  0.55555555555555;
	double w1 =  0.88888888888888;
	double w2 =  0.55555555555555;

LegendrePolynomial(4):
	double x0 = -0.86113631159405; // -sqrt(525 + 70*sqrt(30.0))/35
	double x1 = -0.33998104358486; // -sqrt(525 - 70*sqrt(30.0))/35
	double x2 =  0.33998104358486; //  sqrt(525 - 70*sqrt(30.0))/35
	double x3 =  0.86113631159405; //  sqrt(525 + 70*sqrt(30.0))/35
	double w0 =  0.34785484513745; // (18-sqrt(30.0))/36
	double w1 =  0.65214515486255; // (18+sqrt(30.0))/36
	double w2 =  0.65214515486255; // (18+sqrt(30.0))/36
	double w3 =  0.34785484513745; // (18-sqrt(30.0))/36

LegendrePolynomial(5):
	double x0 = -0.90617984593866; // -sqrt(245 + 14*sqrt(70.0))/21
	double x1 = -0.53846931010568; // -sqrt(245 - 14*sqrt(70.0))/21
	double x2 =  0.00000000000000; //  0
	double x3 =  0.53846931010568; //  sqrt(245 - 14*sqrt(70.0))/21
	double x4 =  0.90617984593866; //  sqrt(245 + 14*sqrt(70.0))/21
	double w0 =  0.236926885056189; // (322 - 13*sqrt(70.0))/900
	double w1 =  0.478628670499366; // (322 + 13*sqrt(70.0))/900
	double w2 =  0.568888888888889; // 128.0/225
	double w3 =  0.478628670499366; // (322 + 13*sqrt(70.0))/900
	double w4 =  0.236926885056189; // (322 - 13*sqrt(70.0))/900
} */
