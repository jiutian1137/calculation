#pragma once
/*{ "clmagic/calculation/algebra/matrix":{
  "Author":"LongJiangnan",
  "Date":"2019-2021",
  "License":"Please identify Author"
} }*/


#include <functional>
#include "utility.h"
#include <iosfwd>
#include <memory>
namespace calculation {
	/* static matrix
	* lattice
	* array
	*/
	template<typename T, size_t M, size_t N, typename ExPo = void>
	class matrixX {
	public:
		using matrix_type = matrixX;
		using scalar_type = T;
		using execution_policy = ExPo;
		
		static constexpr 
		size_t rows() { 
			return M; 
		}
		static constexpr 
		size_t cols() { 
			return N; 
		}
		static constexpr 
		size_t diags() { 
			return (M < N ? M : N); 
		}
		static constexpr 
		size_t rowstep() { 
			return N; 
		}
		static constexpr 
		size_t size() { 
			return M * N; 
		}

		scalar_type _Mydata[rows() * cols()];

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
		
		const scalar_type& at(size_t i, size_t j) const { 
			return _Mydata[i * rowstep() + j]; 
		}
		scalar_type& at(size_t i, size_t j) { 
			return _Mydata[i * rowstep() + j];
		}
		
		/* lattice interface */
		scalar_type operator()(size_t xi, size_t yi) const {
			return at(yi, xi);
		}
		scalar_type operator()(size_t xi, size_t yi, scalar_type xf, scalar_type yf) const {
			auto sample = (*this)(xi,yi);
			if (xf != 0) {
				sample = lerp(sample, (*this)(xi+1,yi), xf);
				if (yf != 0) {
					sample = lerp(sample, lerp((*this)(xi,yi+1),(*this)(xi+1,yi+1),xf), yf);
				}
			} else {
				if (yf != 0) {
					sample = lerp(sample, (*this)(xi,yi+1), yf);
				}
			}

			return sample;
		}
		
		/* array interface */
		const scalar_type& operator[](size_t i) const { 
			return _Mydata[i]; 
		}
		scalar_type& operator[](size_t i) { 
			return _Mydata[i]; 
		}

	public:
		template<typename T2, size_t M2, size_t N2> explicit 
		operator matrixX<T2,M2,N2>() const {
			matrixX<T2,M2,N2> result;
			reshape( MatrixArg4(this->begin(), this->rows(), this->cols(), this->rowstep()),
				MatrixArg4(result.begin(), result.rows(), result.cols(), result.rowstep()) );
			return std::move( result );
		}

		matrix_type operator-() const {
			matrix_type result;
			transform_matrix( 
				/* Left */MatrixArg4(this->begin(), this->rows(), this->cols(), this->rowstep()), 
				/* Result */MatrixArg2(result.begin(), result.rowstep()),
				[](const scalar_type* first, const scalar_type* last, scalar_type* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = -(*first); 
					} 
				});
			return std::move( result );
		}
		
		matrix_type& operator+=(const matrix_type& other) {
			transform_matrix(
				/* Left */MatrixArg4(this->begin(), this->rows(), this->cols(), this->rowstep()), 
				/* Right */MatrixArg2(other.begin(), other.rowstep()),
				/* Result */MatrixArg2(this->begin(), this->rowstep()),
				[](const scalar_type* first1, const scalar_type* last1, const scalar_type* first2, scalar_type* dest) {
					for ( ; first1 != last1; ++first1, ++first2, ++dest) {
						*dest = (*first1) + (*first2);
					}
				});
			return *this;
		}	
		matrix_type operator+(const matrix_type& other) const {
			matrix_type result;
			transform_matrix(
				MatrixArg4(this->begin(), this->rows(), this->cols(), this->rowstep()), 
				MatrixArg2(other.begin(), other.rowstep()),
				MatrixArg2(result.begin(), result.rowstep()),
				[](const scalar_type* first1, const scalar_type* last1, const scalar_type* first2, scalar_type* dest) {
					for ( ; first1 != last1; ++first1, ++first2, ++dest) {
						*dest = (*first1) + (*first2);
					}
				});
			return std::move( result );
		}
		
		matrix_type& operator-=(const matrix_type& other) {
			transform_matrix(
				/* Left */MatrixArg4(this->begin(), this->rows(), this->cols(), this->rowstep()),
				/* Right */MatrixArg2(other.begin(), other.rowstep()),
				/* Result */MatrixArg2(this->begin(), this->rowstep()),
				[](const scalar_type* first1, const scalar_type* last1, const scalar_type* first2, scalar_type* dest) {
					for ( ; first1 != last1; ++first1, ++first2, ++dest) {
						*dest = (*first1) - (*first2);
					}
				});
			return *this;
		}
		matrix_type operator-(const matrix_type& other) const {
			matrix_type result;
			transform_matrix(
				MatrixArg4(this->begin(), this->rows(), this->cols(), this->rowstep()), 
				MatrixArg2(other.begin(), other.rowstep()),
				MatrixArg2(result.begin(), result.rowstep()),
				[](const scalar_type* first1, const scalar_type* last1, const scalar_type* first2, scalar_type* dest) {
					for ( ; first1 != last1; ++first1, ++first2, ++dest) {
						*dest = (*first1) - (*first2);
					}
				});
			return std::move( result );
		}
		
		matrix_type& operator*=(const scalar_type& multiplier) {
			transform_matrix(
				/* Left */MatrixArg4(this->begin(), this->rows(), this->cols(), this->rowstep()),
				/* Result */MatrixArg2(this->begin(), this->rowstep()),
				[&multiplier](const scalar_type* first, const scalar_type* last, scalar_type* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = (*first) * multiplier;
					}
				});
			return *this;
		}
		matrix_type operator*(const scalar_type& multiplier) const {
			matrix_type result;
			transform_matrix(
				/* Left */MatrixArg4(this->begin(), this->rows(), this->cols(), this->rowstep()),
				/* Result */MatrixArg2(result.begin(), result.rowstep()),
				[&multiplier](const scalar_type* first, const scalar_type* last, scalar_type* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = (*first) * multiplier;
					}
				});
			return std::move( result );
		}
		
		matrix_type& operator/=(const scalar_type& divisor) {
			transform_matrix(
				/* Left */MatrixArg4(this->begin(), this->rows(), this->cols(), this->rowstep()),
				/* Result */MatrixArg2(this->begin(), this->rowstep()),
				[&divisor](const scalar_type* first, const scalar_type* last, scalar_type* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = (*first) / divisor;
					}
				});
			return *this;
		}
		matrix_type operator/(const scalar_type& divisor) const {
			matrix_type result;
			transform_matrix(
				/* Left */MatrixArg4(this->begin(), this->rows(), this->cols(), this->rowstep()),
				/* Result */MatrixArg2(result.begin(), result.rowstep()),
				[&divisor](const scalar_type* first, const scalar_type* last, scalar_type* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = (*first) * divisor;
					}
				});
			return std::move( result );
		}
		
		bool operator==(const matrix_type& other) const {
			for (const scalar_type* first2 = other.begin(); const scalar_type & element : *this) {
				if (element != (*first2++)) {
					return false;
				}
			}

			return true;
		}
		
		bool operator!=(const matrix_type& other) const {
			return !(*this == other);
		}
		
		bool operator==(const scalar_type& scalar) const {
			for (const scalar_type& element : *this) {
				if (element != scalar) {
					return false;
				}
			}

			return true;
		}
		
		bool operator!=(const scalar_type& scalar) const {
			return !(*this == scalar);
		}

		template<typename Ty>
		const Ty& as() const { 
			return *reinterpret_cast<const Ty*>(_Mydata); 
		}
		template<typename Ty>
		Ty& as() { 
			return *reinterpret_cast<Ty*>(_Mydata); 
		}

		friend matrix_type operator*(const scalar_type& multiplier, const matrix_type& matrix) {
			return matrix * multiplier;
		}

		friend std::ostream& operator<<(std::ostream& ostr, const matrix_type& matrix) {
			for (size_t i = 0; i != matrix.rows(); ++i) {
				for (size_t j = 0; j != matrix.cols(); ++j) {
					ostr << matrix.at(i,j) << ',';
				}

				ostr << '\n';
			}

			return ostr;
		}
	};

	/* static vector
	* array
	*/
	template<typename T, size_t N, typename ExPo = void>
	using vectorX = matrixX<T, N, 1, ExPo>;

	/* matrix
	* lattice
	* array
	*/
	template<typename T, typename ExPo = void>
	class MatrixX {
	public:
		using matrix_type = MatrixX<T,ExPo>;
		using scalar_type = T;
		using execution_policy = ExPo;

		scalar_type* _Mydata;
		size_t _Mycapacity;
		size_t _Myrows;
		size_t _Mycols;
		size_t _Myrowstep;
		bool _Isref;
	
		size_t capacity() const {
			return _Mycapacity;
		}
		size_t rows() const { 
			return _Myrows; 
		}
		size_t cols() const {
			return _Mycols; 
		}
		size_t diags() const {
			return (_Myrows < _Mycols ? _Myrows : _Mycols); 
		}
		size_t rowstep() const {
			return _Myrowstep; 
		}
		size_t size() const {
			return _Myrows * _Mycols;
		}

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
		
		const scalar_type& at(size_t i, size_t j) const {
			return _Mydata[i * rowstep() + j];
		}
		scalar_type& at(size_t i, size_t j) {
			return _Mydata[i * rowstep() + j];
		}
		
		bool isref() const {
			return _Isref;
		}

		void swap(matrix_type& right) {
			std::swap(_Mydata, right._Mydata);
			std::swap(_Mycapacity, right._Mycapacity);
			std::swap(_Myrows, right._Myrows);
			std::swap(_Mycols, right._Mycols);
			std::swap(_Myrowstep, right._Myrowstep);
			std::swap(_Isref, right._Isref);
		}
		
		void release() {
			if (!_Isref && _Mydata != nullptr) {
				delete[] _Mydata;
			}
			_Mydata = nullptr;
			_Mycapacity = 0;
			_Myrows = 0;
			_Mycols = 0;
			_Myrowstep = 0;
			_Isref = false;
		}
		
		void reallocate(size_t newrows, size_t newcols) {
			if (_Isref) {
				throw std::exception();
			}

			if (_Mydata != nullptr) {
				delete[] _Mydata;
			}
			_Mydata = new scalar_type[newrows * newcols];
			_Mycapacity = newrows * newcols;
			_Myrows = newrows;
			_Mycols = newcols;
			_Myrowstep = _Mycols;
			_Isref = false;
		}

		/* lattice interface */
		scalar_type operator()(size_t xi, size_t yi) const {
			return at(yi, xi);
		}
		scalar_type operator()(size_t xi, size_t yi, scalar_type xf, scalar_type yf) const {
			auto sample = (*this)(xi,yi);
			if (xf != 0) {
				sample = lerp(sample, (*this)(xi+1,yi), xf);
				if (yf != 0) {
					sample = lerp(sample, lerp((*this)(xi,yi+1),(*this)(xi+1,yi+1),xf), yf);
				}
			} else {
				if (yf != 0) {
					sample = lerp(sample, (*this)(xi,yi+1), yf);
				}
			}

			return sample;
		}
		
		/* array interface */
		const scalar_type& operator[](size_t i) const { 
			return _Mydata[i]; 
		}
		scalar_type& operator[](size_t i) { 
			return _Mydata[i]; 
		}

	public:
		MatrixX() : _Mydata(nullptr), _Mycapacity(0), _Myrows(0), _Mycols(0), _Myrowstep(0), _Isref(false) {
			// empty
		}

		MatrixX(MatrixX&& right) noexcept {
			right.swap(*this);
			right.release();
		}

		MatrixX(const MatrixX& right) : MatrixX() {
			this->reallocate(right.rows(), right.cols());
			std::copy(right.begin(), right.end(), this->begin());
		}

		MatrixX(size_t rows, size_t cols) : MatrixX() {
			assert(rows != 0 && cols != 0);
			this->reallocate(rows, cols);
		}

		MatrixX(std::initializer_list<scalar_type> scalar_list) : MatrixX() {
			this->reallocate(scalar_list.size(), 1);
			std::copy(scalar_list.begin(), scalar_list.end(), this->begin());
		}

		MatrixX(size_t rows, size_t cols, std::initializer_list<scalar_type> scalar_list) : MatrixX() {
			assert(rows != 0 && cols != 0);
			assert(rows * cols == scalar_list.size());
			this->reallocate(rows, cols);
			std::copy(scalar_list.begin(), scalar_list.end(), _Mydata);
		}

		~MatrixX() {
			this->release();
		}

		MatrixX& operator=(MatrixX&& right) noexcept {
			right.swap(*this);
			right.release();
			return *this;
		}

		MatrixX& operator=(const MatrixX& right) {
			this->reallocate(right.rows(), right.cols());
			std::copy(right.begin(), right.end(), this->begin());
			return *this;
		}

		MatrixX& operator=(std::initializer_list<scalar_type> scalar_list) {
			if ( isref() ) {
				assert( rows() == scalar_list.size() && cols() == 1 );
				std::copy(scalar_list.begin(), scalar_list.end(), this->begin());
			} else {
				if ( capacity() < scalar_list.size() ) {
					reallocate(scalar_list.size(), 1);
				}

				std::copy(scalar_list.begin(), scalar_list.end(), this->begin());
				_Myrows = scalar_list.size();
				_Mycols = 1;
				_Myrowstep = _Mycols;
			}
			
			return *this;
		}

		void reshape(size_t newrows, size_t newcols) {
			if ( isref() ) {
				throw std::exception();
			}

			size_t newsize = newrows * newcols;
			if (capacity() < newsize) {
				MatrixX newmatrix = MatrixX(newrows, newcols);
				calculation::reshape(MatrixArg4(_Mydata, _Myrows, _Mycols, _Myrowstep), 
					MatrixArg4(newmatrix.data(), newmatrix.rows(), newmatrix.cols(), newmatrix.rowstep()));
				*this = std::move(newmatrix);
			} else {
				calculation::inplace_reshape(MatrixArg4(_Mydata, _Myrows, _Mycols, _Myrowstep),
					newrows, newcols);
				_Myrows = newrows;
				_Mycols = newcols;
				_Myrowstep = _Mycols;
			}
		}

		void fill(const scalar_type& value) {
			std::fill(this->begin(), this->end(), value);
		}

		matrix_type submatrix() {
			abort();
		}

		template<typename T2, typename ExPo2> explicit
		operator MatrixX<T2,ExPo2>() const {
			auto result = MatrixX<T2, ExPo2>(_Myrows, _Mycols);
			calculation::reshape(
				MatrixArg4(_Mydata, _Myrows, _Mycols, _Myrowstep),
				MatrixArg4(result.begin(), result.rows(), result.cols(), result.rowstep()));
			return std::move(result);
		}

		matrix_type& operator+=(const matrix_type& right) {
			assert( _Myrows == right.rows() && _Mycols == right.cols() );
			calculation::transform_matrix(
				MatrixArg4(_Mydata, _Myrows, _Mycols, _Myrowstep),
				MatrixArg2(right.begin(), right.rowstep()),
				MatrixArg2(_Mydata, _Myrowstep),
				[](const scalar_type* first1, const scalar_type* last1, const scalar_type* first2, scalar_type* dest) {
					for ( ; first1 != last1; ++first1, ++first2, ++dest) {
						*dest = (*first1) + (*first2);
					}
				}
			);

			return *this;
		}
		matrix_type operator+(const matrix_type& right) const {
			assert( _Myrows == right.rows() && _Mycols == right.cols() );
			matrix_type result = matrix_type(_Myrows, _Mycols);
			calculation::transform_matrix(
				MatrixArg4(_Mydata, _Myrows, _Mycols, _Myrowstep),
				MatrixArg2(right.begin(), right.rowstep()),
				MatrixArg2(result.begin(), result.rowstep()),
				[](const scalar_type* first1, const scalar_type* last1, const scalar_type* first2, scalar_type* dest) {
					for ( ; first1 != last1; ++first1, ++first2, ++dest) {
						*dest = (*first1) + (*first2);
					}
				}
			);

			return std::move(result);
		}

		matrix_type& operator-=(const matrix_type& right) {
			assert( _Myrows == right.rows() && _Mycols == right.cols() );
			calculation::transform_matrix(
				MatrixArg4(_Mydata, _Myrows, _Mycols, _Myrowstep),
				MatrixArg2(right.begin(), right.rowstep()),
				MatrixArg2(_Mydata, _Myrowstep),
				[](const scalar_type* first1, const scalar_type* last1, const scalar_type* first2, scalar_type* dest) {
					for ( ; first1 != last1; ++first1, ++first2, ++dest) {
						*dest = (*first1) - (*first2);
					}
				}
			);

			return *this;
		}
		matrix_type operator-(const matrix_type& right) const {
			assert( _Myrows == right.rows() && _Mycols == right.cols() );
			matrix_type result = matrix_type(_Myrows, _Mycols);
			calculation::transform_matrix(
				MatrixArg4(_Mydata, _Myrows, _Mycols, _Myrowstep),
				MatrixArg2(right.begin(), right.rowstep()),
				MatrixArg2(result.begin(), result.rowstep()),
				[](const scalar_type* first1, const scalar_type* last1, const scalar_type* first2, scalar_type* dest) {
					for ( ; first1 != last1; ++first1, ++first2, ++dest) {
						*dest = (*first1) - (*first2);
					}
				}
			);

			return std::move(result);
		}
		
		matrix_type& operator*=(const matrix_type& right) {
			assert( _Myrows == _Mycols );
			assert( _Myrows == right.rows() );
			assert( _Mycols == right.cols() );
			calculation::multiply_matrix(
				MatrixArg4(_Mydata, _Myrows, _Mycols, _Myrowstep),
				MatrixArg3(right.begin(), right.cols(), right.rowstep()),
				MatrixArg2(_Mydata, _Myrowstep)
			);

			return *this;
		}
		matrix_type operator*(const matrix_type& right) const {
			assert( _Mycols == right.rows() );
			matrix_type result = matrix_type(_Myrows, right.cols());
			calculation::multiply_matrix(
				MatrixArg4(_Mydata, _Myrows, _Mycols, _Myrowstep),
				MatrixArg3(right.begin(), right.cols(), right.rowstep()),
				MatrixArg2(result.begin(), result.rowstep())
			);

			return std::move(result);
		}
	
		matrix_type& operator*=(const scalar_type& multiplier) {
			calculation::transform_matrix(
				MatrixArg4(_Mydata, _Myrows, _Mycols, _Myrowstep),
				MatrixArg2(_Mydata, _Myrowstep),
				[&multiplier](const scalar_type* first, const scalar_type* last, scalar_type* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = (*first) * multiplier;
					}
				}
			);

			return *this;
		}
		matrix_type operator*(const scalar_type& multiplier) const {
			matrix_type result = matrix_type(_Myrows, _Mycols);
			calculation::transform_matrix(
				MatrixArg4(_Mydata, _Myrows, _Mycols, _Myrowstep),
				MatrixArg2(result.begin(), result.rowstep()),
				[&multiplier](const scalar_type* first, const scalar_type* last, scalar_type* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = (*first) * multiplier;
					}
				}
			);

			return std::move(result);
		}

		matrix_type& operator/=(const scalar_type& divisor) {
			calculation::transform_matrix(
				MatrixArg4(_Mydata, _Myrows, _Mycols, _Myrowstep),
				MatrixArg2(_Mydata, _Myrowstep),
				[&divisor](const scalar_type* first, const scalar_type* last, scalar_type* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = (*first) / divisor;
					}
				}
			);

			return *this;
		}
		matrix_type operator/(const scalar_type& divisor) const {
			matrix_type result = matrix_type(_Myrows, _Mycols);
			calculation::transform_matrix(
				MatrixArg4(_Mydata, _Myrows, _Mycols, _Myrowstep),
				MatrixArg2(result.begin(), result.rowstep()),
				[&divisor](const scalar_type* first, const scalar_type* last, scalar_type* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = (*first) * divisor;
					}
				}
			);

			return std::move(result);
		}

		template<typename Ty>
		const Ty& as() const {
			return *reinterpret_cast<const Ty*>(_Mydata);
		}
		template<typename Ty>
		Ty& as() {
			return *reinterpret_cast<Ty*>(_Mydata);
		}

		friend std::ostream& operator<<(std::ostream& ostr, const matrix_type& matrix) {
			for (size_t i = 0; i != matrix.rows(); ++i) {
				for (size_t j = 0; j != matrix.cols(); ++j) {
					ostr << matrix.at(i,j) << ',';
				}

				ostr << '\n';
			}

			return ostr;
		}
	};

	
#define __calculation_matrix_operator_with_literal(_OP_, _LITERAL_TYPE_)                   \
	template<typename T, size_t M, size_t N, typename ExPo> inline                         \
	matrixX<T,M,N,ExPo> operator##_OP_##(const matrixX<T,M,N,ExPo>& X, _LITERAL_TYPE_ y) { \
		return X _OP_ static_cast<typename matrixX<T,M,N,ExPo>::scalar_type>(y);           \
	}

#define __calculation_matrix_operator_with_literal_commutatibity(_OP_, _LITERAL_TYPE_)     \
	template<typename T, size_t M, size_t N, typename ExPo> inline                         \
	matrixX<T,M,N,ExPo> operator##_OP_##(const matrixX<T,M,N,ExPo>& X, _LITERAL_TYPE_ y) { \
		return X _OP_ static_cast<typename matrixX<T,M,N,ExPo>::scalar_type>(y);           \
	}                                                                                      \
	template<typename T, size_t M, size_t N, typename ExPo> inline                         \
	matrixX<T,M,N,ExPo> operator##_OP_##(_LITERAL_TYPE_ x, const matrixX<T,M,N,ExPo>& Y) { \
		return static_cast<typename matrixX<T,M,N,ExPo>::scalar_type>(x) _OP_ Y;           \
	}

#define __calculation_matrix_lvalueoperator_with_literal(_OP_, _LITERAL_TYPE_)        \
	template<typename T, size_t M, size_t N, typename ExPo> inline                    \
	matrixX<T,M,N,ExPo>& operator##_OP_##(matrixX<T,M,N,ExPo>& X, _LITERAL_TYPE_ y) { \
		return X _OP_ static_cast<typename matrixX<T,M,N,ExPo>::scalar_type>(y);      \
	}

#define __calculation_matrix_comparison_with_literal_commutatibity(_OP_, _LITERAL_TYPE_) \
	template<typename T, size_t M, size_t N, typename ExPo> inline               \
	bool operator##_OP_##(const matrixX<T,M,N,ExPo>& X, _LITERAL_TYPE_ y) {      \
		return X _OP_ static_cast<typename matrixX<T,M,N,ExPo>::scalar_type>(y); \
	}                                                                            \
	template<typename T, size_t M, size_t N, typename ExPo> inline               \
	bool operator##_OP_##(_LITERAL_TYPE_ x, const matrixX<T,M,N,ExPo>& Y) {      \
		return static_cast<typename matrixX<T,M,N,ExPo>::scalar_type>(x) _OP_ Y; \
	}

	__calculation_matrix_operator_with_literal_commutatibity(*, int)
	__calculation_matrix_operator_with_literal_commutatibity(*, float)
	__calculation_matrix_operator_with_literal_commutatibity(*, double)
	__calculation_matrix_operator_with_literal_commutatibity(*, long long)
	__calculation_matrix_operator_with_literal_commutatibity(*, unsigned int)
	__calculation_matrix_operator_with_literal_commutatibity(*, unsigned long long)
	__calculation_matrix_lvalueoperator_with_literal(*=, int)
	__calculation_matrix_lvalueoperator_with_literal(*=, float)
	__calculation_matrix_lvalueoperator_with_literal(*=, double)
	__calculation_matrix_lvalueoperator_with_literal(*=, long long)
	__calculation_matrix_lvalueoperator_with_literal(*=, unsigned int)
	__calculation_matrix_lvalueoperator_with_literal(*=, unsigned long long)
	__calculation_matrix_operator_with_literal(/, int)
	__calculation_matrix_operator_with_literal(/, float)
	__calculation_matrix_operator_with_literal(/, double)
	__calculation_matrix_operator_with_literal(/, long long)
	__calculation_matrix_operator_with_literal(/, unsigned int)
	__calculation_matrix_operator_with_literal(/, unsigned long long)
	__calculation_matrix_lvalueoperator_with_literal(/=, int)
	__calculation_matrix_lvalueoperator_with_literal(/=, float)
	__calculation_matrix_lvalueoperator_with_literal(/=, double)
	__calculation_matrix_lvalueoperator_with_literal(/=, long long)
	__calculation_matrix_lvalueoperator_with_literal(/=, unsigned int)
	__calculation_matrix_lvalueoperator_with_literal(/=, unsigned long long)
	__calculation_matrix_comparison_with_literal_commutatibity(==, int)
	__calculation_matrix_comparison_with_literal_commutatibity(==, float)
	__calculation_matrix_comparison_with_literal_commutatibity(==, double)
	__calculation_matrix_comparison_with_literal_commutatibity(==, long long)
	__calculation_matrix_comparison_with_literal_commutatibity(==, unsigned int)
	__calculation_matrix_comparison_with_literal_commutatibity(==, unsigned long long)
	__calculation_matrix_comparison_with_literal_commutatibity(!=, int)
	__calculation_matrix_comparison_with_literal_commutatibity(!=, float)
	__calculation_matrix_comparison_with_literal_commutatibity(!=, double)
	__calculation_matrix_comparison_with_literal_commutatibity(!=, long long)
	__calculation_matrix_comparison_with_literal_commutatibity(!=, unsigned int)
	__calculation_matrix_comparison_with_literal_commutatibity(!=, unsigned long long)

#undef __calculation_matrix_operator_with_literal
#undef __calculation_matrix_operator_with_literal_commutatibity
#undef __calculation_matrix_comparison_with_literal_commutatibity
	
	template<typename T, size_t M, size_t N, typename ExPo> inline
	matrixX<T,M,N,ExPo> abs(const matrixX<T,M,N,ExPo>& X) {
		if constexpr (M == 0 || N == 0) {
			make_static_array<matrixX<T,M,N,ExPo>>(X, [](const T& x){ return abs(x); });
		} else {
			matrixX<T,M,N,ExPo> result;
			transform_matrix(MatrixArg4(X.begin(), X.rows(), X.cols(), X.rowstep()),
				MatrixArg2(result.begin(), result.rowstep()),
				[](const T* first, const T* last, T* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = abs(*first);
					}
				});
			return std::move( result );
		}
	}

	template<typename T, size_t M, size_t N, typename ExPo> inline
	matrixX<T,M,N,ExPo> floor(const matrixX<T,M,N,ExPo>& X) {
		if constexpr (M == 0 || N == 0) {
			make_static_array<matrixX<T,M,N,ExPo>>(X, [](const T& x){ return floor(x); });
		} else {
			matrixX<T,M,N,ExPo> result;
			transform_matrix(MatrixArg4(X.begin(), X.rows(), X.cols(), X.rowstep()),
				MatrixArg2(result.begin(), result.rowstep()),
				[](const T* first, const T* last, T* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = floor(*first);
					}
				});
			return std::move( result );
		}
	}

	template<typename T, size_t M, size_t N, typename ExPo> inline
	matrixX<T,M,N,ExPo> fract(const matrixX<T,M,N,ExPo>& X) {
		if constexpr (M == 0 || N == 0) {
			make_static_array<matrixX<T,M,N,ExPo>>(X, [](const T& x){ return fract(x); });
		} else {
			matrixX<T,M,N,ExPo> result;
			transform_matrix(MatrixArg4(X.begin(), X.rows(), X.cols(), X.rowstep()),
				MatrixArg2(result.begin(), result.rowstep()),
				[](const T* first, const T* last, T* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = fract(*first);
					}
				});
			return std::move( result );
		}
	}

	template<typename T, size_t M, size_t N, typename ExPo> inline
	matrixX<T,M,N,ExPo> ceil(const matrixX<T,M,N,ExPo>& X) {
		if constexpr (M == 0 || N == 0) {
			make_static_array<matrixX<T,M,N,ExPo>>(X, [](const T& x){ return ceil(x); });
		} else {
			matrixX<T,M,N,ExPo> result;
			transform_matrix(MatrixArg4(X.begin(), X.rows(), X.cols(), X.rowstep()),
				MatrixArg2(result.begin(), result.rowstep()),
				[](const T* first, const T* last, T* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = ceil(*first);
					}
				});
			return std::move( result );
		}
	}

	template<typename T, size_t M, size_t N, typename ExPo> inline
	matrixX<T,M,N,ExPo> round(const matrixX<T,M,N,ExPo>& X) {
		if constexpr (M == 0 || N == 0) {
			make_static_array<matrixX<T,M,N,ExPo>>(X, [](const T& x){ return round(x); });
		} else {
			matrixX<T,M,N,ExPo> result;
			transform_matrix(MatrixArg4(X.begin(), X.rows(), X.cols(), X.rowstep()),
				MatrixArg2(result.begin(), result.rowstep()),
				[](const T* first, const T* last, T* dest) {
					for ( ; first != last; ++first, ++dest) {
						*dest = round(*first);
					}
				});
			return std::move( result );
		}
	}

	template<typename T, size_t M, size_t N, size_t P, typename ExPo> inline
	matrixX<T,M,P,ExPo> operator*(const matrixX<T,M,N,ExPo>& left, const matrixX<T,N,P,ExPo> right) {
		matrixX<T,M,P,ExPo> result;
		multiply_matrix(MatrixArg4(left.begin(), left.rows(), left.cols(), left.rowstep()),
			MatrixArg3(right.begin(), right.cols(), right.rowstep()), 
			MatrixArg2(result.begin(), result.rowstep()));
		return std::move( result );
	}

	template<typename T, size_t N>
	T determinant(const matrixX<T,N,N>& matrix) {
		matrixX<T,N,N> matrix_copy = matrix;

		size_t rowexchange_count;
		if ( matrix_copy.diags() == eliminate(in_marg4(matrix_copy), out_marg2(matrix_copy), &rowexchange_count) ) {
			T det = static_cast<T>(pow( -1, rowexchange_count ));
			for (size_t k = 0; k != matrix_copy.diags(); ++k) {
				det *= matrix_copy.at(k, k);
			}

			return det;
		}
		
		return static_cast<T>(0);
	}

	template<typename T, size_t N>
	matrixX<T,N,N> inverse(const matrixX<T,N,N>& matrix) {
		matrixX<T,N,N> matrix1 = matrix;
		matrixX<T,N,N> matrix2; 
		set_identity( io_marg4(matrix2) );

		size_t rank = 
		eliminate( in_marg4(matrix1), in_marg3(matrix2), 
				   out_marg2(matrix1), out_marg2(matrix2) );
		if ( rank != matrix1.diags() ) { 
			return matrixX<T,N,N>{ 0 }; }
		back_eliminate( in_marg4(matrix1), in_marg3(matrix2), rank, 
						out_marg2(matrix1), out_marg2(matrix2) );
		return matrix2;
	}

	template<typename T, size_t M> inline
	matrixX<T,1,M> transpose(const matrixX<T,M,1>& matrix) {
		return reinterpret_cast<const matrixX<T,1,M>&>(matrix);
	}
	
	template<typename T, size_t N> inline
	matrixX<T,N,1> transpose(const matrixX<T,1,N>& matrix) {
		return reinterpret_cast<const matrixX<T,N,1>&>(matrix);
	}

	/*template<typename T, size_t M, size_t N>
	MatrixX<T,M,N> conjugate_graident_iterate(const SymmetryMatrix_<T,M>& A, const MatrixX<T,M,N>& b, const MatrixX<T,M,N> guess, size_t max_count = 200) {
		auto x = guess;
		auto r = b - A*x;
		auto d = b - A*x;
		for (size_t k = 0; k != max_count; ++k) {
			if ( (r > 0 ? r : -r) <= std::numeric_limits<T>::epsilon() ) {
				break;
			}
			auto alpha = (transpose(r)* r).at(0,0) / (transpose(d) * A * d).at(0,0);
			auto rm1 = r;
			x = x + alpha * d;
			r = r - alpha * A*d;
			double beta = (transpose(r)*r).at(0,0) /(transpose(rm1)*rm1).at(0,0);
			d = r + beta * d;
		}

		return x;
	}*/

	template<typename T, size_t N> inline
	auto dot(const vectorX<T,N>& v1, const vectorX<T,N>& v2) -> typename vectorX<T,N>::scalar_type {
		if constexpr (N == 2) {
			return v1[0]*v2[0] + v1[1]*v2[1];
		} else if constexpr (N == 3) {
			return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
		} else if constexpr (N == 4) {
			return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3];
		} else {
			typename vectorX<T,N>::scalar_type sigma = v1[0] * v2[0];
			for (size_t i = 1; i != N; ++i) {
				sigma += v1[i] * v2[i];
			}

			return sigma;
		}
	}

	template<typename T, size_t N>
	auto norm(const vectorX<T,N>& v1, T L) -> typename vectorX<T, N>::scalar_type {
		typename vectorX<T,N>::scalar_type sigma = pow(v1[0], L);
		for (size_t i = 1; i != N; ++i) {
			sigma += pow(v1[i], L);
		}

		return pow(sigma, static_cast<typename vectorX<T,N>::scalar_type>(1.0/L));
	}

	template<typename T, size_t N> inline
	auto length(const vectorX<T,N>& v1) -> typename vectorX<T,N>::scalar_type {
		return sqrt(dot(v1,v1));
	}
	
	template<typename T, size_t N> inline
	auto distance(const vectorX<T,N>& v1, const vectorX<T,N>& v2) -> typename vectorX<T,N>::scalar_type {
		vectorX<T,N> v3 = v2 - v1;
		return sqrt(dot(v3, v3));
	}

	template<typename T, size_t N> inline
	vectorX<T,N> cross(const vectorX<T,N>& v1, const vectorX<T,N>& v2) {
		/*<idea>
			[ i  j  k  ]
			[ lx ly lz ] = i*1*det(minor(0,0)) + j*-1*det(minor(0,1)) + k*1*det(minor(0,2)), 1.determinat expand
			[ rx ry rz ]
						 = Vector{ det(minor(0,0)), -det(minor(0,1)), det(minor(0,2)) }      2.cast to Vector
		</idea>*/
		return vectorX<T,N>{
			v1[1] * v2[2] - v1[2] * v2[1],
			v1[2] * v2[0] - v1[0] * v2[2],
			v1[0] * v2[1] - v1[1] * v2[0]
		};
	}
	
	template<typename T, size_t N> inline
	vectorX<T,N> cross(const vectorX<T,N>& v1, const vectorX<T,N>& v2, const vectorX<T,N>& v3) {
		/*<idea>
			[   i       j      k     u  ]
			[ v1.x    v1.y   v1.z  v1.w ]
			[ v2.x    v2.y   v2.z  v2.w ] = i*1*det(minor(0,0)) + j*-1*det(minor(0,1)) + k*1*det(minor(0,2)) + u*-1*det(minor(0,3)), 1.determinat expand
			[ v3.x    v3.y   v3.z  v3.w ]
				|      | |    |      |    = Vector{ +(v1.y*detC - v1.z*detE + v1.w*detB),
				+-detA-+-detB-+-detC-+              -(v1.x*detC - v1.z*detF + v1.w*detD),
				|        |    |      |              +(v1.x*detE - v1.y*detF + v1.w*detA),
				+---detD-+----+      |              -(v1.x*detB - v1.y*detD + v1.z*detA) }
				|        |           |
				|   	 +----detE---+
				|                    |
				+-----detF-----------+
		</idea>*/
		T detA = v2[0] * v3[1] - v2[1] * v3[0];
		T detB = v2[1] * v3[2] - v2[2] * v3[1];
		T detC = v2[2] * v3[3] - v2[3] * v3[2];
		T detD = v2[0] * v3[2] - v2[2] * v3[0];
		T detE = v2[1] * v3[3] - v2[3] * v3[1];
		T detF = v2[0] * v3[3] - v2[3] * v3[0];
		return vectorX<T,N>{
			  v1[1]*detC - v1[2]*detE + v1[3]*detB,
			-(v1[0]*detC - v1[2]*detF + v1[3]*detD),
			  v1[0]*detE - v1[1]*detF + v1[3]*detA,
			-(v1[0]*detB - v1[1]*detD + v1[2]*detA) 
		};
	}
	
	template<typename T, size_t N> inline
	vectorX<T,N> normalize(const vectorX<T,N>& v1) {
		T lengthsqr = dot(v1, v1);
		
		// zero vector
		if (lengthsqr == static_cast<T>(0)) {
			return v1;
		}
		// unit vector
		if (lengthsqr == static_cast<T>(1)) {
			return v1;
		}

		return v1 / sqrt(lengthsqr);
	}

	template<typename T>
		using matrix2x2 = matrixX<T, 2, 2>;
	template<typename T>
		using matrix3x3 = matrixX<T, 3, 3>;
	template<typename T>
		using matrix4x4 = matrixX<T, 4, 4>;
	template<typename T>
		using vector2 = vectorX<T, 2>;
	template<typename T>
		using vector3 = vectorX<T, 3>;
	template<typename T>
		using vector4 = vectorX<T, 4>;

	using matrix2x2f = matrixX<float, 2, 2>;
	using matrix3x3f = matrixX<float, 3, 3>;
	using matrix4x4f = matrixX<float, 4, 4>;
	using matrix2x2d = matrixX<double, 2, 2>;
	using matrix3x3d = matrixX<double, 3, 3>;
	using matrix4x4d = matrixX<double, 4, 4>;
	using vector2f = vectorX<float, 2>;
	using vector3f = vectorX<float, 3>;
	using vector4f = vectorX<float, 4>;
	using vector2d = vectorX<double, 2>;
	using vector3d = vectorX<double, 3>;
	using vector4d = vectorX<double, 4>;

	using Matrixf = MatrixX<float>;
	using Matrixd = MatrixX<double>;
}// namespace calculation


#include "../fundamental/array.h"
namespace calculation {
	// { specialization: vector & array  }
	template<typename T, size_t N, typename ExPo>
	class matrixX<T, N, 1, ExPo> {
	public:
		using matrix_type = matrixX;
		using vector_type = matrixX;
		using scalar_type = T;
		using execution_policy = ExPo;

		static constexpr 
		size_t rows() { 
			return N; 
		}
		static constexpr 
		size_t cols() { 
			return 1; 
		}
		static constexpr 
		size_t diags() { 
			return 1; 
		}
		static constexpr 
		size_t rowstep() { 
			return 1; 
		}
		static constexpr 
		size_t size() { 
			return N; 
		}

		scalar_type _Mydata[rows()];

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
		
		const scalar_type& at(size_t i, size_t j) const { 
			assert(j == 0);
			return _Mydata[i]; 
		}
		scalar_type& at(size_t i, size_t j) { 
			assert(j == 0);
			return _Mydata[i]; 
		}
		
		const scalar_type& operator[](size_t i) const { 
			return _Mydata[i]; 
		}
		scalar_type& operator[](size_t i) { 
			return _Mydata[i]; 
		}

	public:
		template<typename T2, size_t M2, size_t N2> explicit
		operator matrixX<T2,M2,N2>() const {
			matrixX<T2,M2,N2> result;
			reshape( MatrixArg4(this->begin(), this->rows(), this->cols(), this->rowstep()),
				MatrixArg4(result.begin(), result.rows(), result.cols(), result.rowstep()) );
			return std::move( result );
		}

		vector_type operator-() const {
			return make_static_array<vector_type>(*this, std::negate<>());
		}
		
		vector_type& operator+=(const vector_type& other) {
			transform_static_array(*this, other, *this, std::plus<>());
			return *this;
		}
		vector_type operator+(const vector_type& other) const {
			return make_static_array<vector_type>(*this, other, std::plus<>());
		}
		
		vector_type& operator-=(const vector_type& other) {
			transform_static_array(*this, other, *this, std::minus<>());
			return *this;
		}
		vector_type operator-(const vector_type& other) const {
			return make_static_array<vector_type>(*this, other, std::minus<>());
		}
		
		vector_type& operator*=(const vector_type& other) {
			transform_static_array(*this, other, (*this), std::multiplies<>());
			return *this;
		}
		vector_type operator*(const vector_type& other) const {
			return make_static_array<vector_type>(*this, other, std::multiplies<>());
		}
		
		vector_type& operator/=(const vector_type& other) {
			transform_static_array(*this, other, (*this), std::divides<>());
			return *this;
		}
		vector_type operator/(const vector_type& other) const {
			return make_static_array<vector_type>(*this, other, std::divides<>());
		}
		
		vector_type& operator%=(const vector_type& other) {
			if constexpr (std::is_same_v<scalar_type, float>) {
				transform_static_array(*this, other, (*this), [](float x, float y){ return _CSTD fmodf(x, y); });
			} else if constexpr (std::is_same_v<scalar_type, double>) {
				transform_static_array(*this, other, (*this), [](double x, double y){ return _CSTD fmod(x, y); });
			} else {
				transform_static_array(*this, other, (*this), std::modulus<>());
			}
			return *this;
		}
		vector_type operator%(const vector_type& other) const {
			// make static array from [*this] and [other] with std::modulus<void>()
			if constexpr (std::is_same_v<scalar_type, float>) {
				return make_static_array<vector_type>(*this, other, [](float x, float y) { return _CSTD fmodf(x, y); });
			} else if constexpr (std::is_same_v<scalar_type, double>) {
				return make_static_array<vector_type>(*this, other, [](double x, double y) { return _CSTD fmod(x, y); });
			} else {
				return make_static_array<vector_type>(*this, other, std::modulus<>());
			}
		}

		vector_type& operator+=(const scalar_type& adder) {
			for (scalar_type& element : *this) {
				element += adder;
			}

			return *this;
		}
		vector_type operator+(const scalar_type& adder) const {
			return make_static_array<vector_type>(*this, [&adder](const scalar_type& x) { 
				return x + adder; });
		}

		vector_type& operator-=(const scalar_type& substractor) {
			for (scalar_type& element : *this) {
				element -= substractor;
			}

			return *this;
		}
		vector_type operator-(const scalar_type& substractor) const {
			return make_static_array<vector_type>(*this, [&substractor](const scalar_type& x) {
				return x - substractor; });
		}

		vector_type& operator*=(const scalar_type& multiplier) {
			for (scalar_type& element : *this) {
				element *= multiplier;
			}

			return *this;
		}
		vector_type operator*(const scalar_type& multiplier) const {
			return make_static_array<vector_type>(*this, [&multiplier](const scalar_type& x) {
				return x * multiplier; });
		}

		vector_type& operator/=(const scalar_type& divisor) {
			for (scalar_type& element : *this) {
				element /= divisor;
			}

			return *this;
		}
		vector_type operator/(const scalar_type& divisor) const {
			return make_static_array<vector_type>(*this, [&divisor](const scalar_type& x) {
				return x / divisor; });
		}

		vector_type& operator%=(const scalar_type& divisor) {
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
		vector_type operator%(const scalar_type& divisor) const {
			if constexpr ( std::is_same_v<scalar_type, float> ) {
				return make_static_array<vector_type>(*this, [&divisor](const scalar_type& x) {
					return _CSTD fmodf(x, divisor); });
			} else if constexpr ( std::is_same_v<scalar_type, double> ) {
				return make_static_array<vector_type>(*this, [&divisor](const scalar_type& x) {
					return _CSTD fmod(x, divisor); });
			} else {
				return make_static_array<vector_type>(*this, [&divisor](const scalar_type& x) {
					return x % divisor; });
			}
		}

		bool operator==(const matrix_type& other) const {
			for (const scalar_type* first2 = other.begin(); const scalar_type& element : *this) {
				if (element != (*first2++)) {
					return false;
				}
			}

			return true;
		}
		
		bool operator!=(const matrix_type& other) const {
			return !(*this == other);
		}
		
		bool operator==(const scalar_type& scalar) const {
			for (const scalar_type& element : *this) {
				if (element != scalar) {
					return false;
				}
			}

			return true;
		}
		
		bool operator!=(const scalar_type& scalar) const {
			return !(*this == scalar);
		}

		template<typename Ty>
		const Ty& as() const { 
			return *reinterpret_cast<const Ty*>(_Mydata); 
		}
		template<typename Ty>
		Ty& as() { 
			return *reinterpret_cast<Ty*>(_Mydata); 
		}

		friend vector_type operator+(const scalar_type& adder, const vector_type& vector) {
			return vector + adder;
		}
		friend vector_type operator-(const scalar_type& substractor, const vector_type& vector) {
			return vector - substractor;
		}
		friend vector_type operator*(const scalar_type& multiplier, const vector_type& vector) {
			return vector * multiplier;
		}
	};

#define __calculation_vector_operator_with_literal(_OP_, _LITERAL_TYPE_)               \
	template<typename T, size_t N, typename ExPo> inline                               \
	vectorX<T,N,ExPo> operator##_OP_##(const vectorX<T,N,ExPo>& X, _LITERAL_TYPE_ y) { \
		return X _OP_ static_cast<typename vectorX<T,N,ExPo>::scalar_type>(y);         \
	}

#define __calculation_vector_operator_with_literal_commutatibity(_OP_, _LITERAL_TYPE_) \
	template<typename T, size_t N, typename ExPo> inline                               \
	vectorX<T,N,ExPo> operator##_OP_##(const vectorX<T,N,ExPo>& X, _LITERAL_TYPE_ y) { \
		return X _OP_ static_cast<typename vectorX<T,N,ExPo>::scalar_type>(y);         \
	}                                                                                  \
	template<typename T, size_t N, typename ExPo> inline                               \
	vectorX<T,N,ExPo> operator##_OP_##(_LITERAL_TYPE_ x, const vectorX<T,N,ExPo>& Y) { \
		return static_cast<typename vectorX<T,N,ExPo>::scalar_type>(x) _OP_ Y;         \
	}

#define __calculation_vector_lvalueoperator_with_literal(_OP_, _LITERAL_TYPE_)    \
	template<typename T, size_t N, typename ExPo> inline                          \
	vectorX<T,N,ExPo>& operator##_OP_##(vectorX<T,N,ExPo>& X, _LITERAL_TYPE_ y) { \
		return X _OP_ static_cast<typename vectorX<T,N,ExPo>::scalar_type>(y);    \
	}

	__calculation_vector_operator_with_literal_commutatibity(+, int)
	__calculation_vector_operator_with_literal_commutatibity(+, float)
	__calculation_vector_operator_with_literal_commutatibity(+, double)
	__calculation_vector_operator_with_literal_commutatibity(+, long long)
	__calculation_vector_operator_with_literal_commutatibity(+, unsigned int)
	__calculation_vector_operator_with_literal_commutatibity(+, unsigned long long)
	__calculation_vector_lvalueoperator_with_literal(+=, int)
	__calculation_vector_lvalueoperator_with_literal(+=, float)
	__calculation_vector_lvalueoperator_with_literal(+=, double)
	__calculation_vector_lvalueoperator_with_literal(+=, long long)
	__calculation_vector_lvalueoperator_with_literal(+=, unsigned int)
	__calculation_vector_lvalueoperator_with_literal(+=, unsigned long long)
	__calculation_vector_operator_with_literal(-, int)
	__calculation_vector_operator_with_literal(-, float)
	__calculation_vector_operator_with_literal(-, double)
	__calculation_vector_operator_with_literal(-, long long)
	__calculation_vector_operator_with_literal(-, unsigned int)
	__calculation_vector_operator_with_literal(-, unsigned long long)
	__calculation_vector_lvalueoperator_with_literal(-=, int)
	__calculation_vector_lvalueoperator_with_literal(-=, float)
	__calculation_vector_lvalueoperator_with_literal(-=, double)
	__calculation_vector_lvalueoperator_with_literal(-=, long long)
	__calculation_vector_lvalueoperator_with_literal(-=, unsigned int)
	__calculation_vector_lvalueoperator_with_literal(-=, unsigned long long)

#undef __calculation_vector_operator_with_literal
#undef __calculation_vector_operator_with_literal_commutatibity
#undef __calculation_vector_lvalueoperator_with_literal
}// vector array_operator


// undeterminant
//namespace calculation {
//	// { find matrix pivot position }
//	template<typename _MatTy>
//	struct const_pivot_iterator {
//		using matrix_type = _MatTy;
//		using Scalar = matrix_scalar_t<_MatTy>;
//
//		static std::array<size_t,2> _First_major(const matrix_type& _Matrix) {
//			// search a not equal than 0.0+-_Threshould in every colume from {0,0}
//			std::array<size_t,2> _Where = { 0, 0 };
//			for ( ; _Where[1] != _Matrix.cols(); ++_Where[1], _Where[0] = 0) {
//				for ( ; _Where[0] != _Matrix.rows(); ++_Where[0]) {// Has major in the colume ?
//					if ( abs(_Matrix.at(_Where[0],_Where[1])) > std::numeric_limits<Scalar>::epsilon() ) {
//						return _Where;
//					}
//				}
//			}
//
//			return std::array<size_t,2>{ _Matrix.rows(), _Matrix.cols() };
//		}
//
//		static std::array<size_t,2> _Last_major(const matrix_type& _Matrix) {
//			// inverse_search a not equal than 0.0+-_Threshould in every row from {_Matrix.rows(), _Matrix.cols()}
//			std::array<size_t,2> _Where = { _Matrix.rows(), _Matrix.cols() };
//			do {
//				_Where[0] -= 1;
//				_Where[1] = 0;
//				for ( ; _Where[1] != _Matrix.cols(); ++_Where[1]) {// Has major in the row ?
//					if (abs(_Matrix.at(_Where[0],_Where[1])) > std::numeric_limits<Scalar>::epsilon() ) {
//						return _Where;
//					}
//				}
//			} while (_Where[0] != 0);
//
//			return std::array<size_t,2>{ static_cast<size_t>(-1) , static_cast<size_t>(-1) };
//		}
//
//		static std::array<size_t,2> _Next_major(const matrix_type& _Matrix, std::array<size_t,2> _Pos) {// _Next major-pos, _Pos must be valid major-pos
//			assert( _Pos[0] < _Matrix.rows() && _Pos[1] < _Matrix.cols() );
//
//			std::array<size_t,2> _Where = { _Pos[0] + 1, _Pos[1] + 1 };
//			for ( ; _Where[1] != _Matrix.cols(); ++_Where[1], _Where[0] = _Pos[0] + 1) {
//				for ( ; _Where[0] != _Matrix.rows(); ++_Where[0]) {// [ _Pos[1]+1, rows )
//					if ( abs(_Matrix.at(_Where[0],_Where[1])) > std::numeric_limits<Scalar>::epsilon() ) {
//						return _Where;
//					}
//				}
//			}
//
//			return std::array<size_t,2>{ _Matrix.rows(), _Matrix.cols() };
//		}
//
//		static std::array<size_t,2> _Prev_major(const matrix_type& _Matrix, std::array<size_t,2> _Pos) {// _Pos must be valid major-pos
//			assert( _Pos[0] < _Matrix.rows() && _Pos[1] < _Matrix.cols() );
//
//			if (_Pos[0] == _Pos[1]) {// Rank(_Matrix) == _Matrix.diags()
//				return std::array<size_t, 2>{ _Pos[0] - 1, _Pos[1] - 1 };
//			} else {
//				std::array<size_t,2> _Where = { _Pos[0] - 1, _Pos[1] - 1 };
//				do {
//					_Where[0] -= 1;
//					_Where[1] = 0;
//					for ( ; _Where[1] != _Pos[1]; ++_Where[1]) {// Has major in the row ?
//						if ( abs(_Matrix.at(_Where[0],_Where[1])) > std::numeric_limits<Scalar>::epsilon() ) {
//							return _Where;
//						}
//					}
//				} while (_Where[0] != 0);
//				
//				return std::array<size_t,2>{ static_cast<size_t>(-1), static_cast<size_t>(-1) };
//			}
//		}
//	
//		explicit const_pivot_iterator(const matrix_type& _Matrix) : matrix_ptr(&_Matrix), position(_First_major(_Matrix)) { 
//			// Check major position
//			for (size_t i = 0; i != position[0]; ++i) {
//				for (size_t j = position[1] + 1; j != matrix_ptr->cols(); ++j) {
//					if ( abs(matrix_ptr->at(i, j)) > std::numeric_limits<Scalar>::epsilon() ) {
//						throw std::exception("calculation::pivot_iterator::no-swapable");
//					}
//				}
//			}
//		}
//
//		const_pivot_iterator(const matrix_type& _Matrix, std::array<size_t,2> _Mjpos) : matrix_ptr(&_Matrix), position(_Mjpos) {
//			/*if (_Mjpos[0] == _Matrix.rows() && _Mjpos[1] == _Matrix.cols()) {
//		
//			} else {
//				if ( _Mjpos[0] > _Matrix.rows() || _Mjpos[1] > _Matrix.cols() ) {
//					throw std::exception("clmagic::major_iterator::invalid-major-pos");
//				}
//				if ( abs(_Matrix.at(_Mjpos[0], _Mjpos[1])) < std::numeric_limits<Scalar>::epsilon() ) {
//					throw std::exception("clmagic::major_iterator::invalid-major-pos");
//				}
//			}*/
//		}
//
//		const Scalar& operator*() const {
//			return matrix_ptr->at(position[0], position[1]);
//		}
//
//		const_pivot_iterator& operator++() {
//			// 1. Get _Next major position
//			std::array<size_t, 2> _Next = _Next_major(*matrix_ptr, position);
//			// 2. Check _Next major position
//			for (size_t i = position[0] + 1; i != _Next[0]; ++i) {
//				for (size_t j = position[1] + 1; j != matrix_ptr->cols(); ++j) {
//					if (abs(matrix_ptr->at(i, j)) > std::numeric_limits<Scalar>::epsilon()) {
//						throw std::exception("calculation::pivot_iterator::no-swapable");
//					}
//				}
//			}
//
//			position = _Next;
//			return (*this);
//		}
//		const_pivot_iterator& operator--() {
//			position = _Prev_major(*matrix_ptr, position);
//			return (*this);
//		}
//		const_pivot_iterator operator--(int) {
//			const_pivot_iterator _Tmp = *this;
//			--(*this);
//			return _Tmp;
//		}
//		const_pivot_iterator operator++(int) {
//			const_pivot_iterator _Tmp = *this;
//			++(*this);
//			return _Tmp;
//		}
//
//		bool operator==(const const_pivot_iterator& _Right) const {
//			return (matrix_ptr == _Right.matrix_ptr) && (position[0] == _Right.position[0]) && (position[1] == _Right.position[1]);
//		}
//		bool operator!=(const const_pivot_iterator& _Right) const {
//			return !(*this == _Right);
//		}
//		bool operator<(const const_pivot_iterator& _Right) const {
//			assert(matrix_ptr == _Right.matrix_ptr);
//			return (position[0] + position[1]) < (_Right.position[0] + _Right.position[1]);
//		}
//		bool operator>(const const_pivot_iterator& _Right) const {
//			return (_Right < *this);
//		}
//		bool operator<=(const const_pivot_iterator& _Right) const {
//			return !(*this > _Right);
//		}
//		bool operator>=(const const_pivot_iterator& _Right) const {
//			return !(*this < _Right);
//		}
//
//		const Scalar* ptr() const {
//			return matrix_ptr->data() + position[0] * matrix_ptr->cols() + position[1];
//		}
//		size_t size() const {
//			return matrix_ptr->cols() - position[1];
//		}
//
//		void seek_to_first() {
//			this->position = _First_major(*matrix_ptr);
//		}
//		void seek_to_last() {
//			this->position = _Last_major(*matrix_ptr);
//		}
//
//		const matrix_type* matrix_ptr;
//		std::array<size_t,2> position;
//	};
//
//	// { find matrix pivot position }
//	template<typename _MatTy>
//	struct pivot_iterator : public const_pivot_iterator<_MatTy> {
//		using _Mybase     = const_pivot_iterator<_MatTy>;
//		using matrix_type = _MatTy;
//		using Scalar = matrix_scalar_t<_MatTy>;
//
//		explicit pivot_iterator(matrix_type& _Matrix) : _Mybase(_Matrix) {}
//
//		pivot_iterator(matrix_type& _Matrix, std::array<size_t,2> _Mjpos) : _Mybase(_Matrix, _Mjpos) {}
//
//		Scalar& operator*() {
//			return const_cast<Scalar&>(_Mybase::operator*());
//		}
//
//		pivot_iterator& operator++() {
//			std::array<size_t,2> _Next = _Mybase::_Next_major(*_Mybase::matrix_ptr, _Mybase::position);
//			if (_Next[0] != _Mybase::position[0] + 1 && /*safe_check*/_Next[0] < _Mybase::matrix_ptr->rows()) {
//				matrix_type& _Matrix = const_cast<matrix_type&>(*_Mybase::matrix_ptr);
//				matrix_row_swap(_Matrix, _Next[0], _Matrix, _Mybase::position[0] + 1);
//				_Next[0] = _Mybase::position[0] + 1;
//			}
//
//			_Mybase::position = _Next;
//			return *this;
//		}
//
//		pivot_iterator& operator--() {
//			_Mybase::operator--();
//			return *this;
//		}
//
//		pivot_iterator operator--(int) {
//			pivot_iterator _Tmp = *this;
//			--(*this);
//			return _Tmp;
//		}
//
//		pivot_iterator operator++(int) {
//			pivot_iterator _Tmp = *this;
//			++(*this);
//			return _Tmp;
//		}
//	};
//}// namespace undeterminant

/*<undetermined>
	template<typename genScalar> inline
	void minor(matrix_slice<const genScalar> _Source, matrix_slice<genScalar> _Dest, size_t i, size_t j) {
		size_t _Rseek = 0, _Cseek;
		for (size_t _Rfirst = 0; _Rfirst < _Dest.rows(); ++_Rfirst) {
			if (_Rfirst != i) {
				_Cseek = 0;
				for (size_t _Cfirst = 0; _Cfirst < _Dest.cols(); ++_Cfirst) {
					if (_Cfirst != j) {
						_Dest.at(_Rseek, _Cseek) = _Source.at(_Rfirst, _Cfirst);
						++_Cseek;
					}
				}
				++_Rseek;
			}
		}
	}
</undetermined>*/

#ifdef _INCLUDED_MM2
//namespace calculation {
//	template<size_t _Rows, size_t _Cols>
//	using m128matrix = matrix<float, _Rows, _Cols, __m128>;
//
//	using m128matrix4x4 = matrix<float, 4, 4, __m128>;
//}
//
//inline calculation::m128matrix<4,4> operator+(const calculation::m128matrix<4,4>& matrix1, const calculation::m128matrix<4,4>& matrix2) {
//    calculation::m128matrix<4, 4> result;
//    const __m128* lhs = reinterpret_cast<const __m128*>(matrix1.data());
//    const __m128* rhs = reinterpret_cast<const __m128*>(matrix2.data());
//    __m128*       dst = reinterpret_cast<__m128*>(result.data());
//
//    *dst++ = _mm_add_ps(*lhs++, *rhs++);
//    *dst++ = _mm_add_ps(*lhs++, *rhs++);
//    *dst++ = _mm_add_ps(*lhs++, *rhs++);
//    *dst = _mm_add_ps(*lhs, *rhs);
//
//    return std::move(result);
//}
//
//inline calculation::m128matrix<4,4> operator-(const calculation::m128matrix<4,4>& matrix1, const calculation::m128matrix<4,4>& matrix2) {
//    calculation::m128matrix<4, 4> result;
//    const __m128* lhs = reinterpret_cast<const __m128*>(matrix1.data());
//    const __m128* rhs = reinterpret_cast<const __m128*>(matrix2.data());
//    __m128*       dst = reinterpret_cast<__m128*>(result.data());
//
//    *dst++ = _mm_sub_ps(*lhs++, *rhs++);
//    *dst++ = _mm_sub_ps(*lhs++, *rhs++);
//    *dst++ = _mm_sub_ps(*lhs++, *rhs++);
//    *dst = _mm_sub_ps(*lhs, *rhs);
//
//    return std::move(result);
//}
//
//inline calculation::m128matrix<4,4> operator*(const calculation::m128matrix<4,4>& matrix, const float scalar) {
//    calculation::m128matrix<4, 4> result;
//    const __m128* lhs = reinterpret_cast<const __m128*>(matrix.data());
//    const __m128  rhs = _mm_set1_ps(scalar);
//    __m128*       dst = reinterpret_cast<__m128*>(result.data());
//
//    *dst++ = _mm_sub_ps(*lhs++, rhs);
//    *dst++ = _mm_sub_ps(*lhs++, rhs);
//    *dst++ = _mm_sub_ps(*lhs++, rhs);
//    *dst = _mm_sub_ps(*lhs, rhs);
//
//    return std::move(result);
//}
//
//inline calculation::m128matrix<4,4> operator*(const calculation::m128matrix<4,4>& matrix1, const calculation::m128matrix<4,4>& matrix2) {
//    calculation::m128matrix<4,4> result;
//    const float*  lhs = matrix1.data();
//    const __m128* rhs = reinterpret_cast<const __m128*>(matrix2.data());
//    __m128*       dst = reinterpret_cast<__m128*>(result.data());
//
//    for (size_t i = 0; i != 4; ++i) {
//        *dst = _mm_mul_ps(_mm_set1_ps(*lhs), *rhs); ++lhs; ++rhs;
//        *dst = _mm_add_ps(*dst, _mm_mul_ps(_mm_set1_ps(*lhs), *rhs)); ++lhs; ++rhs;
//        *dst = _mm_add_ps(*dst, _mm_mul_ps(_mm_set1_ps(*lhs), *rhs)); ++lhs; ++rhs;
//        *dst = _mm_add_ps(*dst, _mm_mul_ps(_mm_set1_ps(*lhs), *rhs)); ++lhs; rhs = reinterpret_cast<const __m128*>(matrix2.data()); ++dst;
//    }
//
//    return std::move(result);
//}
//
//inline calculation::m128matrix<4,4> transpose(const calculation::m128matrix<4,4>& matrix) {
//    calculation::m128matrix<4,4> result;
//
//    __m128 _Tmp3, _Tmp2, _Tmp1, _Tmp0;
//    _Tmp0 = _mm_shuffle_ps(matrix.at<__m128>(0,0), matrix.at<__m128>(1,0), 0x44);
//    _Tmp2 = _mm_shuffle_ps(matrix.at<__m128>(0,0), matrix.at<__m128>(1,0), 0xEE);
//    _Tmp1 = _mm_shuffle_ps(matrix.at<__m128>(2,0), matrix.at<__m128>(3,0), 0x44);
//    _Tmp3 = _mm_shuffle_ps(matrix.at<__m128>(2,0), matrix.at<__m128>(3,0), 0xEE);
//
//    result.at<__m128>(0,0) = _mm_shuffle_ps(_Tmp0, _Tmp1, 0x88);
//    result.at<__m128>(1,0) = _mm_shuffle_ps(_Tmp0, _Tmp1, 0xDD);
//    result.at<__m128>(2,0) = _mm_shuffle_ps(_Tmp2, _Tmp3, 0x88);
//    result.at<__m128>(3,0) = _mm_shuffle_ps(_Tmp2, _Tmp3, 0xDD);
//
//    return std::move(result);
//}
//
//inline calculation::m128matrix<4,4> inverse(const calculation::m128matrix<4,4>& matrix) {
//    auto result = inverse(reinterpret_cast<const calculation::matrix<float,4,4>&>(matrix));
//    return reinterpret_cast<const calculation::m128matrix<4,4>&>(result);
//}

#endif // _INCLUDED_MM2
