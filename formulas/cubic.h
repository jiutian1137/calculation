/*{ "clmagic/calculation/formulas/cubic":{
  "Mathematician": ["Cardano","Bombelli","François Viète"],
  "Reference": [
    // Cardano
    "https://math.vanderbilt.edu/schectex/courses/cubic/",
    // Viète
    "http://www.sosmath.com/algebra/factor/fac11/fac11.html",
    "https://mathworld.wolfram.com/CubicFormula.html"
  ]
  "Author": "LongJiangnan",
  "Date": "2021",
  "License": "Please identify Mathematician",
} }*/
#pragma once


#include <assert.h>
#include <math.h>
//#include <complex>
namespace calculation {
	// Polynomial: a*pow(x,3) + b*pow(x,2) + c*pow(x,1) + d
	template<typename Real>
	struct Cubic {
		Real a, b, c, d;
		Real operator()(Real x) const {
			return ((a*x + b)*x + c)*x + d;
		}
	};

	// d/dx * ( a*pow(x,3) + b*pow(x,2) + c*pow(x,1) + d )
	template<typename Real> inline
	Real differentiate(Cubic<Real> f, Real x) {
		return 3 * f.a * x + 2 * f.b * x + f.c * x;
	}

	// Solve a*pow(z,3) + b*pow(z,2) + c*pow(z,1) + d = 0 respect 'z'
	template<typename Real>
	size_t solve_zeros(Cubic<Real> polynomial, Real* roots) {
		assert(roots != nullptr);
		auto [a, b, c, d] = polynomial;

		/* Eliminate "a":{
			pow(z,3) + b/a*pow(z,2) + c/a*pow(z,1) + d/a = 0
		} */
		d /= a;
		c /= a;
		b /= a;
		//a /= a;

		/* 
		Define:{
			        b
			z = x - -
			        3
		},
		Eliminate Any*pow(z,2):{
			               3*c - b*b                  9*b*c - 27*d - 2*b*b*b
			pow(x,3) + 3 * --------- * pow(x,1) - 2 * ---------------------- = 0
			                   9                               54
		},
		Simplify:{
			pow(x,3) + p*pow(x,1) - q = 0
		}*/
		Real p = (3*c - b*b) / 3;
		Real q = (9*b*c - 27*d - 2*b*b*b) / 27;

		/* 
		Vieta's subtitution:{ 
					 p
			x = w - ---
					3*w
		},
		Subtitute and simplify:{
						  pow(p,3)
			pow(w,3) - ----------- - q = 0
					   27*pow(w,3)
		},
		Multiply by "pow(w,3)":{
					     pow(p,3)
			pow(w,3*2) - -------- - q*pow(w,3) = 0
						    27
		},
		Solve the quadratic:{
			           1              1       1
			pow(w,3) = - * q +- sqrt( -*q*q + --*p*p*p ) 
			           2              4       27
		} */
		Real R = q / 2;
		Real Q = p / 3;
		Real QQQ = Q * Q * Q;
		Real discrim = R * R + QQQ;
		if (discrim < 0) {
			Real pi = 3.141592653589793;
			Real theta = acos(R / sqrt(-QQQ));
			roots[0] = 2 * sqrt(-Q) * cos((theta - 2*pi) / 3) - b / 3;
			roots[1] = 2 * sqrt(-Q) * cos(theta / 3) - b / 3;
			roots[2] = 2 * sqrt(-Q) * cos((theta + 2*pi) / 3) - b / 3;
			return 3;
		}

		// Subtitute
		Real sqrt_discrim = sqrt(discrim);
		Real w = cbrt(R + sqrt_discrim);
		Real x = w - p/(3*w);
		Real z = x - b/3;
		roots[0] = z;
		return 1;
	}

	// Gerolamo-Cardano's method
}// namespace cubic