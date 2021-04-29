#pragma once
/*{ "clmagic/calculation/formulas/bernstein":{
  "Author":"LongJiangnan",
  "Date":"2019-2021",
  "License":"Please identify Mathmatician"
} }*/


#include <assert.h>
#include "binomial.h"
namespace calculation {
	using BernsteinPolynomialOrder = size_t;
	
	struct BernsteinBasis {
		BernsteinPolynomialOrder n;
		size_t i;

		template<typename Real>
		Real operator()(Real t) const {
			return binomial_coefficient(n,i) * pow(t, i)* pow(1 - t, n - i);
		}
	};

	struct BernsteinPolynomial {
		BernsteinPolynomialOrder n;

		BernsteinBasis operator[](size_t i) const {
			return BernsteinBasis{n,i};
		}

		// sum(0,n, Cn(i) * pow(t,i) * pow(1-t,n-i))
		template<typename Real>
		Real operator()(Real t) const {
			auto Cn = binomial_coefficient(n);
			Real result = 0;
			for (size_t i = 0; i <= n; ++i) {
				result += Cn[i] * pow(t, i) * pow(1-t, n-i);
			}

			return std::move(result);
		}

		// sum(0,n, f(i) * basis)
		template<typename Real, typename Fn>
		auto operator()(Real t, Fn f) const {
			auto Cn = binomial_coefficient(n);

			Real basis  = Cn[0]/* * pow(t, 0)*/ * pow(1-t, n/*-0*/);
			auto result = f(0) * basis;

			for (size_t i = 1; i <= n; ++i) {
				basis  = Cn[i] * pow(t, i) * pow(1-t, n-i);
				result += f(i) * basis;
			}

			return std::move(result);
		}

		// sum(0,n, x[i] * basis)
		template<typename Real, typename Iter>
		auto operator()(Real t, Iter x, Iter x_last) const {
			assert(std::distance(x, x_last) == n+1);
			auto Cn = binomial_coefficient(n);

			Real basis  = Cn[0] * pow(1-t, n);
			auto result = (*x++) * basis;
			
			for (size_t i = 1; i <= n; ++i) {
				basis  = Cn[i] * pow(t, i) * pow(1-t, n-i);
				result += (*x++) * basis;
			}

			return std::move(result);
		}
	};

	// d/dt * B_n_i
	template<typename Real> inline
	Real differentiate(BernsteinBasis B_n_i, Real t) {
		const auto [n, i] = B_n_i;
		return n * (BernsteinBasis{n-1,i-1}(t) - BernsteinBasis{n-1,i}(t));
	}
}