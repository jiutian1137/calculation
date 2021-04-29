#pragma once
/*{ "clmagic/calculation/formulas/quadratic": {
  "Author": "LongJiangnan",
  "Date": "2021",
  "License": "Free"
} }*/


#include <math.h>// sqrt(..)
namespace calculation {
	/*{ Polynomial: a*pow(x,2) + b*pow(x,1) + c*pow(x,0) ,
		circle_equation ,
		sphere_equation ,
	}*/
	template<typename Real>
	struct Quadratic {
		Real a, b, c;
		Real operator()(Real x) const {
			return (a*x + b)*x + c;
		}
	};

	template<typename Real> inline
	Real differentiate(Quadratic<Real> f, Real x) {
		// d/dx * ( a*pow(x,2) + b*pow(x,1) + c*pow(x,0) )
		return 2 * f.a * x + f.b;
	}
	
	template<typename Real>
	size_t solve_zeros(Quadratic<Real> polynomial, decltype(__nullptr)) {
		const auto& [a, b, c] = polynomial;
		Real discrim = b * b - a * c * 4;
		return discrim >= 0 ? 2 : 0;
	}

	template<typename Real>
	size_t solve_zeros(Quadratic<Real> polynomial, Real* roots) {
		// notation
		const auto& [a, b, c] = polynomial;
		// check has roots
		Real discrim = b * b - a * c * 4;
		if ( discrim >= 0 ) {
			// switch numeric situation
			if ( b > 0 ) {
				Real negb_sub_sqrtdiscrim = -b - sqrt(discrim);
				roots[0] = negb_sub_sqrtdiscrim / (a * 2);
				roots[1] = (c * 2) / negb_sub_sqrtdiscrim;
			} else {
				Real negb_add_sqrtdiscrim = -b + sqrt(discrim);
				roots[0] = (c * 2) / negb_add_sqrtdiscrim;
				roots[1] = negb_add_sqrtdiscrim / (a * 2);
			}
		}

		return discrim >= 0 ? 2 : 0;
	}
}// namespace quadratic

/* {
"Principle": {
	when b < 0
	root1 = [-b + sqrt(b*b-4*a*c)] / (2*a)
	root0 = [-b - sqrt(b*b-4*a*c)] / (2*a)

				[-b - sqrt(b*b-4*a*c)] * [-b + sqrt(b*b-4*a*c)]
			= --------------------------------------------------------     : multiple [-b + sqrt(b*b-4*a*c)]
									2*a * [-b + sqrt(b*b-4*a*c)]

				(-b)*(-b) - (b*b - 4*a*c)	             4*a*c
			= ---------------------------- = ----------------------------  : eliminate
			2*a * [-b + sqrt(b*b-4*a*c)]   2*a * [-b + sqrt(b*b-4*a*c)]

			= 2*c/[-b + sqrt(b*b-4*a*c)]
},
"Reference": { Book, Numerical Analysis, Timothy-Sauer }
} */