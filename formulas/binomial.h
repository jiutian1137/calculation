#pragma once
/*{ "clmagic/calculation/formulas/binomial":{
  "Author":"LongJiangnan",
  "Date":"2021",
  "License":"Please identify Author and Mathmatician"
} }*/


#include <vector>
namespace calculation {
	// positive order
	using BinomialOrder = size_t;

	// { faster than recursion more 100 times }
	template<typename Real = size_t>
	Real binomial_coefficient(BinomialOrder N, size_t i) {
		// C(n,0) = C(n,n) = 1
		if (i == 0 || i == N) {
			return 1;
		}

		// C(n,1) = C(n,n-1) = n
		if (i == 1 || i + 1 == N) {
			return N;
		}

		// symmetry
		// ...

		std::vector<Real> cache;
		cache.resize(i+1, 1);
		for (BinomialOrder n = 2; n != N; ++n) {
			// clip right{ ignore_terms }
			size_t k = std::min(size_t(n), i/*cache.size()-1*/);
			// clip right{ C(n,n) }
			if (k == size_t(n)) {
				--k;
			}
			// clip left{ ignore_terms }
			size_t final_k = i/*cache.size()-1*/ - std::min(N-n, i);
			// clip left{ C(n,0) }
			size_t last_k = final_k == 0 ? 0 : final_k - 1;
						
			for ( ; k != last_k; --k) {
				cache[k] += cache[k - 1];
			}
		}

		return cache[i] + cache[i - 1];
	}

	// { binomial coefficients of order N }
	template<typename Real = size_t>
	std::vector<Real> binomial_coefficient(BinomialOrder N) {
		std::vector<Real> cache;
		cache.resize(N+1, 1);
		if (N >= 2) {
			for (BinomialOrder n = 2; n <= N; ++n) {
				// clip right{ C(n,n) }
				size_t k = n - 1;
				// clip left{ C(n,0) }
				size_t last_k = 0;
			
				for ( ; k != last_k; --k) {
					cache[k] += cache[k - 1];
				}
			}
		}

		return std::move(cache);
	}

	// { set lower_triangular of MatrixCoef to be binomial_coefficent with n order, factorial(n) / (factorial(n-1) * factorial(i)) }
	template<typename Matrix>
	void set_positive_binomial_coefficient_matrix( Matrix& C, uint32_t n ) {
		assert( n >= 0 && n+1 <= C.diags() );

		C.at(0,0) = 1;
		if (n == 0) {
			return;
		}

		for (size_t i = 1; i <= n; ++i) {
			C.at(i,0) = 1;
			C.at(i,i) = 1;
			for (size_t k = 1; k != i; ++k) {
				C.at(i,k) = C.at(i-1,k) + C.at(i-1,k-1);
			}
		}
	}

	// { C.at(-i,k) = pow(-1,k) * C.at(i+k-1,k) }
	template<typename Matrix, typename InMatrix>
	void set_negative_binomial_coefficient_matrix( Matrix& nC, int n, const InMatrix& C ) {
		/*for (size_t i = 1; i <= abs(n); ++i) {
			for (size_t k = 1; k != i; ++k) {
				nC.at(i,k) = pow(-1,k) * C.at(i+k-1, k);
			}
		}*/
	}
}