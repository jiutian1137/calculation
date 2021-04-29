#pragma once
#include <cmath>
#include <numeric>
#include <numbers>

namespace calculation {
	using ChebyshevPolynomialOrder = int;

	/* T(n) = cos(n*acos(x)) */
	struct ChebyshevPolynomialFirstKindTag {  
		ChebyshevPolynomialOrder order;
		ChebyshevPolynomialFirstKindTag(ChebyshevPolynomialOrder L) : order(L) {}
	};

	template<typename Real>
	Real evaluate(ChebyshevPolynomialFirstKindTag T, Real x) {
		if (T.order == 0) {
			return static_cast<Real>(1);
		}

		Real Tprev = static_cast<Real>(1);
		Real Tcurrent = x;
		for (ChebyshevPolynomialOrder i = 1; i < T.order; ++i) {
			Real Tnext = Tcurrent * ( 2 * x ) - Tprev;
			Tprev = Tcurrent;
			Tcurrent = Tnext;
		}
		return Tcurrent;
	}

	template<typename Real>
	size_t solve_zeros(ChebyshevPolynomialFirstKindTag T, Real* roots) {
		// { n*acos(x) = odd*(pi/2), x = cos(odd*(pi/2n)) }
		// if odd = 2n, x = cos(pi), next repeat.
		if (roots) {
			const Real pi = std::numbers::pi_v<Real>;
			const uint32_t n = static_cast<uint32_t>(T.order);

			for (uint32_t odd = 1; odd <= 2*n; odd += 2) {
				(*roots++) = cos( (pi * odd) / (n * 2) );
			}
		}

		return static_cast<size_t>(T.order);
	}

	// ...

}// namespace calculation