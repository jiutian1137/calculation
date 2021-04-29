/*{ "clmagic/calculation/algebra/utility":{ 
  "Author":"LongJiangnan",
  "Date":"2020-2021",
  "License":"Please identify Author"
} }*/
#pragma once


#include <concepts>
#include <type_traits>
#include <iterator>
namespace calculation {
	template<typename Type>
	concept vectorconcept = requires(Type v) {
		v.size();
		v.begin();
		v.end();
		v[(size_t)(2)];
	};

	template<typename Type>
	concept matrixconcept = requires(Type m) {
		m.rows();
		m.cols();
		m.diags();
		m.rowstep();

		m.begin();
		m.end();
		m.at((size_t)(4), (size_t)(4));
	};

	template<matrixconcept Matrix>
	using matrix_scalar_t = std::remove_cv_t<std::remove_reference_t< decltype( Matrix().at((size_t)(4),(size_t)(4)) ) >>;
}


#include <cassert>
#include <tuple>
#include <vector>// diagmatrix_fixed_point_iterate
#include <numeric>
#include <algorithm>
namespace calculation { 
	/*template<typename Fn1, typename Fn2, typename Fn3>
	void submatrix_for(Fn1 break_if_condition, Fn2 next_submatrix, Fn3 body) {
		while (true) {
			body();
			if (break_if_condition()) {
				break;
			}
			next_submatrix();
		}
	}*/
	
	template<typename Iter>
	struct MatrixArg4 {
		Iter first;
		size_t rows;
		size_t cols;
		size_t rowstep;
		
		MatrixArg4() = delete;

		MatrixArg4(Iter first, size_t rows, size_t cols, size_t rowstep)
			: first(first), rows(rows), cols(cols), rowstep(rowstep) {}
	};

	template<typename Iter>
	struct MatrixArg3 {
		Iter first;
		size_t cols;
		size_t rowstep;

		MatrixArg3() = delete;

		MatrixArg3(Iter first, size_t cols, size_t rowstep)
			: first(first), cols(cols), rowstep(rowstep) {}

		MatrixArg3(MatrixArg4<Iter> args)
			: first(args.first), cols(args.cols), rowstep(args.rowstep) {}
	};

	template<typename Iter>
	struct MatrixArg2 {
		Iter first;
		size_t rowstep;

		MatrixArg2() = delete;

		MatrixArg2(Iter first, size_t rowstep) : first(first), rowstep(rowstep) {}

		MatrixArg2(MatrixArg4<Iter> args) : first(args.first), rowstep(args.rowstep) {}
	};


	template<matrixconcept InMatrix> inline
	auto in_marg4(const InMatrix& matrix) -> MatrixArg4<decltype(matrix.begin())> {
		return MatrixArg4<decltype(matrix.begin())>(
			matrix.begin(), 
			matrix.rows(), 
			matrix.cols(), 
			matrix.rowstep()
			);
	}
	
	template<matrixconcept InMatrix> inline
	auto out_marg4(InMatrix& matrix) -> MatrixArg4<decltype(matrix.begin())> {
		return MatrixArg4<decltype(matrix.begin())>(
			matrix.begin(), 
			matrix.rows(), 
			matrix.cols(), 
			matrix.rowstep()
			);
	}

	template<matrixconcept InMatrix> inline
	auto io_marg4(InMatrix& matrix) -> MatrixArg4<decltype(matrix.begin())> {
		return MatrixArg4<decltype(matrix.begin())>(
			matrix.begin(), 
			matrix.rows(), 
			matrix.cols(), 
			matrix.rowstep()
			);
	}

	template<matrixconcept InMatrix> inline
	auto in_marg3(const InMatrix& matrix) -> MatrixArg3<decltype(matrix.begin())> {
		return MatrixArg3<decltype(matrix.begin())>(
			matrix.begin(), 
			matrix.cols(), 
			matrix.rowstep()
			);
	}
	
	template<matrixconcept InoutMatrix> inline
	auto out_marg3(InoutMatrix& matrix) -> MatrixArg3<decltype(matrix.begin())> {
		return MatrixArg3<decltype(matrix.begin())>(
			matrix.begin(), 
			matrix.cols(), 
			matrix.rowstep()
			);
	}

	template<matrixconcept InoutMatrix> inline
	auto io_marg3(InoutMatrix& matrix) -> MatrixArg3<decltype(matrix.begin())> {
		return MatrixArg3<decltype(matrix.begin())>(
			matrix.begin(), 
			matrix.cols(), 
			matrix.rowstep()
			);
	}

	template<matrixconcept InMatrix> inline
	auto in_marg2(const InMatrix& matrix) -> MatrixArg2<decltype(matrix.begin())> {
		return MatrixArg2<decltype(matrix.begin())>(
			matrix.begin(), 
			matrix.rowstep()
			);
	}

	template<matrixconcept InoutMatrix> inline
	auto out_marg2(InoutMatrix& matrix) -> MatrixArg2<decltype(matrix.begin())> {
		return MatrixArg2<decltype(matrix.begin())>(
			matrix.begin(), 
			matrix.rowstep()
			);
	}

	template<matrixconcept InoutMatrix> inline
	auto io_marg2(InoutMatrix& matrix) -> MatrixArg2<decltype(matrix.begin())> {
		return MatrixArg2<decltype(matrix.begin())>(
			matrix.begin(), 
			matrix.rowstep()
			);
	}

	template<typename Iter>
	struct VectorArrayArg4 {
		Iter first;
		size_t count;
		size_t vsize;
		size_t stride;
		VectorArrayArg4(Iter first, size_t count, size_t vsize, size_t stride)
			: first(first), count(count), vsize(vsize), stride(stride) {}
	};



	template<typename Ty>
	struct Ones {
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = intptr_t;
		using value_type = Ty;
		using reference = Ty&;
		using pointer = Ty*;

		Ty value;

		Ones() = default;
		Ones(const Ty& val) : value(val) {}

		const Ty& operator[](size_t) const {
			return value;
		}
		Ty& operator[](size_t) {
			return value;
		}

		const Ty& operator*() const {
			return value;
		}
		Ty& operator*() {
			return value;
		}

		Ones& operator++() {
			return *this;
		}
		Ones& operator--() {
			return *this;
		}
		
		Ones& operator++(int) {
			return *this;
		}
		Ones& operator--(int) {
			return *this;
		}

		Ones& operator+=(difference_type diff) {
			return *this;
		}
		Ones& operator-=(difference_type diff) {
			return *this;
		}

		Ones operator+(difference_type diff) const {
			return *this;
		}
		Ones operator-(difference_type diff) const {
			return *this;
		}

		// no comparator
	};

	// set psource[i] into destination[i]
	template<typename InputIterator, typename OutputIterator>
	void set_matrix(InputIterator psource, MatrixArg4<OutputIterator> destination) {
		using Scalar = typename std::iterator_traits<OutputIterator>::value_type;
		while ( true ) {
			// set first_row
			for (size_t i = 0; i != destination.cols; ++i, 
				++destination.first, ++psource) {
				*destination.first = static_cast<Scalar>( *psource );
			}

			// break if final_row
			if ( destination.rows == 1 ) {
				break;
			}

			// remove destination[0]
			std::advance(destination.first, 
				destination.rowstep - destination.cols);
			--destination.rows;
		}
	}

	// set psource[i] into destination[i]
	template<typename InputIterator, matrixconcept Matrix> inline
	void set_matrix(InputIterator psource, Matrix& destination) {
		set_matrix(psource, MatrixArg4(destination.begin(), 
			destination.rows(), destination.cols(), destination.rowstep()));
	}

	// set psource[i] into destination.row(k).element(i)
	template<typename InputIterator, typename OutputIterator>
	void set_row(InputIterator psource, MatrixArg4<OutputIterator> destination, size_t k) {
		// [k]th_row to first_row
		std::advance(destination.first, k*destination.rowstep);
		destination.rows -= k;

		// copy the row
		using Scalar = typename std::iterator_traits<OutputIterator>::value_type;
		auto dest      = destination.first;
		auto dest_last = std::next(destination.first, destination.cols);
		for ( ; dest != dest_last; ++dest, ++psource) {
			*dest = static_cast<Scalar>( *psource );
		}
	}

	// set psource[i] into destination.row(k).element(i)
	template<typename InputIterator, matrixconcept Matrix> inline
	void set_row(InputIterator psource, Matrix& destination, size_t k) {
		set_row(psource, MatrixArg4(destination.begin(),
			destination.rows(), destination.cols(), destination.rowstep()),
			k);
	}

	// set psource[i] into destination.column(k).element(i)
	template<typename InputIterator, typename OutputIterator>
	void set_column(InputIterator psource, MatrixArg4<OutputIterator> destination, size_t k) {
		// [k]th_column to first_column
		std::advance(destination.first, k);
		destination.cols -= k;

		// copy the colume
		using Scalar = typename std::iterator_traits<OutputIterator>::value_type;
		while (true) {
			*destination.first = static_cast<Scalar>( *psource );
			
			if (destination.rows == 1) {
				break;
			}

			std::advance(destination.first, destination.rowstep);
			--destination.rows;
			++psource;
		}
	}
	
	// set psource[i] into destination.column(k).element(i)
	template<typename InputIterator, matrixconcept Matrix> inline
	void set_column(InputIterator psource, Matrix& destination, size_t k) {
		set_column(psource, MatrixArg4(destination.begin(),
			destination.rows(), destination.cols(), destination.rowstep()),
			k);
	}

	// set psource[i] into destination.diagonal(k).element(i)
	template<typename InputIterator, typename OutputIterator>
	void set_diagonal(InputIterator psource, MatrixArg4<OutputIterator> destination, int k) {
		// [k]th_diagonal to main_diagonal
		if (k > 0) {
			std::advance(destination.first, k);
			destination.cols -= k;
		} else if(k < 0) {
			std::advance(destination.first, (-k) * destination.rowstep);
			destination.rows -= (-k);
		}

		// copy the diagonal
		using Scalar = typename std::iterator_traits<OutputIterator>::value_type;
		size_t destination_diags = std::min(destination.rows, destination.cols);
		while (true) {
			*destination.first = static_cast<Scalar>( *psource++ );

			if (destination_diags == 1) {
				break;
			}

			std::advance(destination.first, destination.rowstep + 1);
			--destination_diags;
			++psource;
		}
	}
	
	// set psource[i] into destination.diagonal(k).element(i)
	template<typename InputIterator, matrixconcept Matrix> inline
	void set_diagonal(InputIterator psource, Matrix& destination, int k) {
		set_diagonal(psource, MatrixArg4(destination.begin(),
			destination.rows(), destination.cols(), destination.rowstep()),
			k);
	}

	// set '1' into destination.diagonal(0).element(i)
	template<typename OutputIterator> inline
	void set_matrix_identity(MatrixArg4<OutputIterator> destination) {
		set_matrix(Ones(0), destination);
		set_diagonal(Ones(1), destination, 0);
	}
	
	// set '1' into destination.diagonal(0).element(i)
	template<matrixconcept Matrix> inline
	void set_matrix_identity(Matrix& destination) {
		set_matrix_identity(MatrixArg4(destination.begin(),
			destination.rows(), destination.cols(), destination.rowstep()));
	}

	

	template<typename InputIterator, typename OutputIterator, typename Fn>
	void transform_matrix(MatrixArg4<InputIterator> source, MatrixArg2<OutputIterator> destination, Fn ranges_f) {
		while ( true ) {
			// do transform in the row
			ranges_f(
				source.first, std::next(source.first, source.cols),
				destination.first
			);

			// break if final_row transformed
			if ( source.rows == 1 ) {
				break;
			}

			// remove row[0]
			std::advance(source.first, source.rowstep);
			std::advance(destination.first, destination.rowstep);
			--source.rows;
		}
	}

	template<typename InputMatrix, typename OutputMatrix, typename Fn> inline
	void transform_matrix(const InputMatrix& source, OutputMatrix& destination, Fn ranges_f) {
		transform_matrix(MatrixArg4(source.begin(), source.rows(), source.cols(), source.rowstep()),
			MatrixArg2(destination.begin(), destination.rowstep()),
			ranges_f);
	}

	template<typename InputIterator1, typename InputIterator2, typename OutputIterator, typename Fn>
	void transform_matrix(MatrixArg4<InputIterator1> matrix1, MatrixArg2<InputIterator2> matrix2, MatrixArg2<OutputIterator> result, Fn ranges_f) {
		while ( true ) {
			// do transform in the row
			ranges_f(
				matrix1.first, std::next(matrix1.first, matrix1.cols),
				matrix2.first, 
				result.first
			);

			// break if final_row transformed
			if ( matrix1.rows == 1 ) {
				break;
			}

			// remove row[0]
			std::advance(matrix1.first, matrix1.rowstep);
			std::advance(matrix2.first, matrix2.rowstep);
			std::advance(result.first, result.rowstep);
			--matrix1.rows;
		}
	}
	 
	template<matrixconcept InputMatrix1, matrixconcept InputMatrix2, matrixconcept OutputMatrix, typename Fn> inline
	void transform_matrix(const InputMatrix1& matrix1, const InputMatrix2& matrix2, OutputMatrix& result, Fn ranges_f) {
		transform_matrix(MatrixArg4(matrix1.begin(), matrix1.rows(), matrix1.cols(), matrix1.rowstep()),
			MatrixArg2(matrix2.begin(), matrix2.rowstep()),
			MatrixArg2(result.begin(), result.rowstep()),
			ranges_f);
	}



	struct _LinearCombination {
		// { linear combination, sum(vectors[i] * weights[i]) = out_vector }
		template<typename InputIterator1, typename InputIterator2, typename OutputIterator>
		void operator()(VectorArrayArg4<InputIterator1> vector_array, InputIterator2 weight_array, OutputIterator out_vector) const {
			assert( vector_array.count != 0 );
			OutputIterator out_vector_last = std::next(out_vector, vector_array.vsize);
		
			// out_vector = vector_array[0] * weight_array[0]
			OutputIterator out_vector_i = out_vector;
			InputIterator1 vector_i     = vector_array.first;
			const auto& weight = *weight_array;
			for ( ; out_vector_i != out_vector_last; ++out_vector_i, ++vector_i) {
				*out_vector_i = (*vector_i) * weight;
			}
			if ( vector_array.count == 1 ) {
				return;
			}

			do {
				std::advance(vector_array.first, vector_array.stride);
				vector_array.count -= 1;
				std::advance(weight_array, 1);

				// out_vector += vector_array[i] * weight_array[i], i below [1,...)
				out_vector_i = out_vector;
				vector_i     = vector_array.first;
				const auto& weight = *weight_array;
				for ( ; out_vector_i != out_vector_last; ++out_vector_i, ++vector_i) {
					*out_vector_i += (*vector_i) * weight;
				}
			} while ( vector_array.count != 1 );
		}
	};

	/* {
		"Function": multiply_matrix(MatrixMxN, MatrixNxP, MatrixMxP, LinearCombOperator),
		"Application": { "MatrixTransform", "MatrixLinearCombination" , "Span" },
		"Relation": multiply_matrix2x2(A,B,ref(C)),
		"Principle": {
			left_colume_vectors and right_colume_weights;
			| V1x| V2x| V3x |   | w1 ... |
			| V1y| V2y| V3y | * | w2 ... | = | V1*w1 + V2*w2 + V3*w3 |
			| V1z| V2z| V3z |   | w3 ... |
			equal
			right_row_vectors and left_row_weights, and row_memory_continuous
					       | V1x V1y ... | 
			| w1 w2 w3 | * | V2x V2y ... | = | w1*V1 + w2*V2 + w3*V3 |
						   | V3x V3y ... | 
		}
	} */
	template<typename InputIterator1, typename InputIterator2, typename OutputIterator, typename Fn>
	void multiply_matrix(MatrixArg4<InputIterator1> matrix1, MatrixArg3<InputIterator2> matrix2, MatrixArg2<OutputIterator> result, Fn __lincomb) {
		while (true) {
			// linear_combination_of( row_vectors[i] * weights[i] )
			auto row_vectors = VectorArrayArg4<InputIterator2>(matrix2.first, matrix1.cols, matrix2.cols, matrix2.rowstep);
			InputIterator1 weights    = matrix1.first;
			OutputIterator out_vector = result.first;
			__lincomb(row_vectors, weights, out_vector);

			// break if final_combination
			if ( matrix1.rows == 1 ) {
				break;
			}
			
			// next combination
			std::advance(matrix1.first, matrix1.rowstep);
			std::advance(result.first, result.rowstep);
			--matrix1.rows;
		}
	}

	// matrix1(MxN) * matrix2(NxP) = result(MxP)
	template<typename InputIterator1, typename InputIterator2, typename OutputIterator> inline
	void multiply_matrix(MatrixArg4<InputIterator1> matrix1, MatrixArg3<InputIterator2> matrix2, MatrixArg2<OutputIterator> result) {
		multiply_matrix( matrix1, matrix2, result, _LinearCombination() );
	}

	// Strassen matrix multiplication, Matrix2x2 * Matrix2x2 = Matrix2x2
	template<typename FwdIt1, typename FwdIt2, typename OutputIterator>
	void multiply_matrix2x2(MatrixArg2<FwdIt1> A, MatrixArg2<FwdIt2> B, MatrixArg2<OutputIterator> C) {
		const auto& A11 = *A.first++;
		const auto& A12 = *A.first; std::advance(A.first, A.rowstep - 1);
		const auto& A21 = *A.first++;
		const auto& A22 = *A.first;
		const auto& B11 = *B.first++;
		const auto& B12 = *B.first; std::advance(B.first, B.rowstep - 1);
		const auto& B21 = *B.first++;
		const auto& B22 = *B.first;
		auto& C11 = *C.first++;
		auto& C12 = *C.first; std::advance(C.first, C.rowstep - 1);
		auto& C21 = *C.first++;
		auto& C22 = *C.first;

		auto P1 = A11 * (B12 - B22);
		auto P2 = (A11 + A12) * B22;
		auto P3 = (A21 + A22) * B11;
		auto P4 = A22 * (B21 - B11);
		auto P5 = (A11 + A22) * (B11 + B22);
		auto P6 = (A12 - A22) * (B21 + B22);
		auto P7 = (A11 - A21) * (B11 + B12);
		C11 = P5 + P4 - P2 + P6;
		C12 = P1 + P2;
		C21 = P3 + P4;
		C22 = P5 + P1 - P3 - P7;
	}



	// source(MxN) to destination(AxB)
	template<typename InputIterator, typename OutputIterator>
	void reshape(MatrixArg4<InputIterator> source, MatrixArg4<OutputIterator> destination) {
		using OutputTraits = std::iterator_traits<OutputIterator>;
		using Scalar = typename OutputTraits::value_type;

		size_t copied_cols = std::min(source.cols, destination.cols);
		while ( true ) {
			// copy and fillzero in first_row
			size_t j = 0;
			for ( ; j != copied_cols; ++j) {
				*destination.first++ = static_cast<Scalar>(*source.first++);
			}
			for ( ; j != destination.cols; ++j) {
				*destination.first++ = static_cast<Scalar>(0);
			}

			// break if (source or destination is final_row)  
			if ( source.rows == 1 || destination.rows == 1 ) {
				// fill zeros to destination and break
				if ( destination.rows != 1 ) {
					std::advance(destination.first, destination.rowstep - destination.cols);
					--destination.rows;
					set_matrix(Ones(0), destination);
				}

				break;
			}

			// next submatrix(rows-1,cols)
			std::advance(source.first, source.rowstep - copied_cols);
			--source.rows;
			std::advance(destination.first, destination.rowstep - destination.cols);
			--destination.rows;
		}
	}
	
	// source(MxN) to destination(AxB)
	template<matrixconcept InputMatrix, matrixconcept OutputMatrix> inline
	void reshape(const InputMatrix& source, OutputMatrix destination) {
		reshape(MatrixArg4(source.begin(), source.rows(), source.cols(), source.rowstep()),
			MatrixArg4(destination.begin(), destination.rows(), destination.cols(), destination.rowstep()));
	}

	// source(MxN) to source(AxB)
	template<typename BidirIterator>
	void inplace_reshape(MatrixArg4<BidirIterator> matrix, size_t newrows, size_t newcols) {
		if ( newcols <= matrix.cols ) {
			reshape(matrix, MatrixArg4(matrix.first, newrows, newcols, newcols));
		} else {
			using Scalar = typename std::iterator_traits<BidirIterator>::value_type;
			
			size_t newrowstep = newcols;
			BidirIterator source_back = std::next(matrix.first,
				(std::min(matrix.rows, newrows) - 1) * matrix.rowstep + matrix.cols - 1);
			BidirIterator destination_back = std::next(matrix.first, 
				(std::min(matrix.rows,newrows) - 1) * newrowstep + newcols - 1);
			size_t copied_cols = std::min(matrix.cols, newcols);
			
			if ( newrows > matrix.rows ) {
				set_matrix( 
					Ones(0), 
					MatrixArg4(
						std::next(destination_back,newrowstep-std::max(matrix.cols,newcols)+1),
						newrows - matrix.rows,
						newcols,
						newrowstep
					));
				newrows = matrix.rows;
			}

			while ( true ) {
				size_t j = newcols;
				for ( ; j != copied_cols; --j) {
					*destination_back-- = static_cast<Scalar>(0);
				}
				for ( ; j != 0; --j) {
					*destination_back-- = static_cast<Scalar>(*source_back--);
				}

				if (newrows == 1) {
					break;
				}

				// next submatrix(rows-1,cols)
				std::advance(source_back, -ptrdiff_t(matrix.rowstep - copied_cols));
				--matrix.rows;
				std::advance(destination_back, -ptrdiff_t(newrowstep - newcols));
				--newrows;
			}
		}
	}

	// source(MxN) to source(AxB)
	template<matrixconcept Matrix> inline
	void inplace_reshape(Matrix& matrix, size_t newrows, size_t newcols, size_t newrowstep) {
		inplace_reshape(MatrixArg4(matrix.begin(), matrix.rows(), matrix.cols(), matrix.rowstep()),
			newrows, newcols, newrowstep);
	}

	// source(MxN) to destination(NxM)
	template<typename InputIterator, typename OutputIterator>
	void transpose(MatrixArg4<InputIterator> source, MatrixArg2<OutputIterator> destination) {
		while ( true ) {
			// copy row(source) and colume(destination)
			OutputIterator transposed_colume = destination.first;
			for (size_t k = 0, final_k = source.cols-1; true; ) {
				*transposed_colume = *source.first;
				if ( k == final_k ) {
					break;
				}
				std::advance(source.first, 1);
				std::advance(transposed_colume, destination.rowstep);
				++k;
			}

			// break if only row has been tranposed
			if (source.rows == 1) {
				break;
			}
			
			// next submatrix, source(rows-1,cols) destination(rows,cols-1)
			std::advance(source.first, source.rowstep - (source.cols-1));
			--source.rows;
			std::advance(destination.first, 1);
		}
	}
	
	// source(MxN) to destination(NxM)
	template<matrixconcept InputMatrix, matrixconcept OutputMatrix> inline
	void transpose(const InputMatrix& source, OutputMatrix destination) {
		transpose( MatrixArg4(source.begin(), source.rows(), source.cols(), source.rowstep()),
			MatrixArg2(destination.begin(), destination.rowstep()) );
	}

	/* {
		"Function": transpose_inplace(first, rows, cols, rowstep),
		"Principle": {
			replace cols with rowstep,
			next_index = (index * rows) % size + index / cols,
			next_index = (index * rows) % (size - 1),

			create swap_cycles,
			apply swap_cycles
		}
	} */
	template<typename FwdIt>
	void inplace_transpose(MatrixArg4<FwdIt> matrix) {
		std::vector<std::vector<size_t>> swap_cycles;
		size_t size = matrix.rows * matrix.rowstep;
		size_t counter = 0;
		
		// create swap_cycles
		for (size_t i = 0; i < matrix.rows && counter != size-1; ++i) {
			for (size_t j = 0; j < matrix.cols && counter != size-1; ++j) {
				size_t index = i * matrix.rowstep + j;

				// break if element swaped
				if ( !swap_cycles.empty() ) {
					auto first = swap_cycles.begin();
					auto last  = swap_cycles.end();
					bool swaped = false;
					do {
						swaped |= std::find(first->begin(), first->end(), index) != first->end();
					} while (++first != last && !swaped);

					if (swaped) {
						continue;
					}
				}

				// create a new swap_cycle
				swap_cycles.push_back(std::vector<size_t>());
				size_t next_index = index;
				do {
					swap_cycles.back().push_back(next_index);
					next_index = (next_index * matrix.rows) % size
								+ next_index / matrix.rowstep;
				} while (next_index != index);

				// counter swaped elements
				counter += swap_cycles.back().size();
			}
		}

		// apply swap_cycles
		std::for_each(swap_cycles.begin(), swap_cycles.end(), 
			[matrix](std::vector<size_t>& swap_cycle) {
				auto _First = swap_cycle.begin();
				auto _Last  = swap_cycle.end();
				++_First;
				for ( ; _First != _Last; ++_First) {
					std::swap(*std::next(matrix.first, swap_cycle.front()),
								*std::next(matrix.first, *_First));
				}
			});

		//for (size_t i = 1; i != size_m1; ++i, ++first) {
		//	if (flags[i]) {
		//		continue;
		//	}
		//	size_t j = (i * rows) % size_m1;// i below [1,size-1]
		//	// size_t j = (i*rows)%size + i/cols;// i below [0,size]
		//	std::swap( *first, *std::next(ori_first, j) );
		//	flags[i] = flags[j] = true;
		//}
	}

	template<matrixconcept Matrix> inline
	void inplace_transpose(Matrix& matrix) {
		inplace_transpose(MatrixArg4(matrix.begin(),
			matrix.rows(), matrix.cols(), matrix.rowstep()));
	}


	/*_INFORMATION({
		"Function": eliminate_matrix(...),
		"Relation": [ bareiss_eliminate(...), eliminate_augumented_matrix(...) ],
		"Application": this is part of solving equation, another part is back-subtitute,
		"Principle": { 
			"Standard": {
				for (i to rows-1) {
					find pivot;
					rowexchange first and finded, and counter;
					eliminate lower rows;
				}
			},
			"Method1": {{
				saved origin_rows;
				zero rowexchange_count;
				while (true) {
					find pivot;
					break if not finded pivot;
					break if submatrix only one row;
					remove zero_columns;
					rowexchange first and finded, and counter;
					eliminate lower rows;
					break if only column has been processed;
					remove row[0] and remove column[0];
				}
			}, difference: we only added 2 steps and some conditions }
		}
	})*/

	/*_INFORMATION({
		"Principle": {
			| 1           |
			| 0 1 0       | this is full_rank_matrix, but is singular_matrix;
			|       1     |

			| 1           |
			| 0 1         | this is non_singular_matrix;
			| 0 0 1       |

			(Singular)ignore lower rows;
			goto back pivot;
			while (true) {
				find(test) pivot;
				(Singular)accelerate;
				scale pivot row;
				break if the only row has been scaled;
				eliminate upper rows;
				remove row[0], prev pivot;
			}
		}
	})*/

	struct _Find_pivot {
		template<typename FwdIt>
		std::tuple<bool,FwdIt> operator()(MatrixArg4<FwdIt> matrix) const {
			while (true) {
				// find pivot for colume
				FwdIt seek = matrix.first;
				for (size_t i = 0, final_i = matrix.rows-1; true; ) {
					if ( *seek != 0 ) {
						return { true, seek };
					}
					if ( i == final_i ) {
						break;
					}
					std::advance(seek, matrix.rowstep);
					++i;
				}

				// break if the only column has been processed
				if (matrix.cols == 1) {
					break;
				}

				// remove column[0]
				std::advance(matrix.first, 1);
				--matrix.cols;
			}

			return { false, matrix.first };
		}
	};

	struct _Find_maxabs_pivot {
		template<typename FwdIt>
		std::tuple<bool,FwdIt> operator()(MatrixArg4<FwdIt> matrix) const {
			while ( true ) {
				// find pivot for colume
				FwdIt seek = matrix.first;
				FwdIt pivot = matrix.first;
				for (size_t i = 0, final_i = matrix.rows-1; true; ) {
					if( abs(*seek) > abs(*pivot) ) {
						pivot = seek;
					}
					if ( i == final_i ) {
						break;
					}

					std::advance(seek, matrix.rowstep);
					++i;
				}
				if ( *pivot != 0 ) {
					return { true, pivot };
				}

				// break if the only column has been processed
				if ( matrix.cols == 1 ) {
					break;
				}

				// remove column[0]
				std::advance(matrix.first, 1);
				--matrix.cols;
			}

			return { false, matrix.first };
		}
	};

	struct _Swap_ranges {
		template<typename Iter1, typename Iter2>
		void operator()(Iter1 first_row, Iter1 first_row_last, Iter2 finded_row) const {
			std::swap_ranges(first_row, first_row_last, finded_row);
		}
	};

	struct _Fma_ranges {
		/* [first, last) * multiplier + [dest, ...) = [dest, ...) */
		template<typename InputIterator, typename Scl, typename OutputIterator>
		void operator()(InputIterator first, InputIterator last, const Scl& multiplier, OutputIterator dest) const {
			for ( ; first != last; ++first, ++dest) {
				*dest += (*first) * multiplier;
			}
		}
	};

	struct _Scale_ranges {
		/* [first, last) * multiplier = [dest, ...) */
		template<typename InputIterator, typename Ty, typename OutputIterator>
		void operator()(InputIterator first, InputIterator last, const Ty& multiplier, OutputIterator dest) const {
			for ( ; first != last; ++first, ++dest) {
				*dest = (*first) * multiplier;
			}
		}
	};

	struct _Dscale_ranges {
		/* [first, last) / divisor = [dest, ...) */
		template<typename InputIterator, typename Ty, typename OutputIterator>
		void operator()(InputIterator first, InputIterator last, const Ty& divisor, OutputIterator dest) const {
			for ( ; first != last; ++first, ++dest) {
				*dest = (*first) / divisor;
			}
		}
	};
	
	template<typename FwdIt, typename Fn1, typename Fn2, typename Fn3>
	size_t inplace_eliminate(MatrixArg4<FwdIt> matrix, size_t* rowexchange_count, 
		Fn1 find_pivot, Fn2 swap_ranges, Fn3 fma_ranges) {
		assert( matrix.rows != 0 && matrix.cols != 0 );

		// saved origin_rows
		size_t origin_rows = matrix.rows;

		// zero rowexchange_count
		if (rowexchange_count) {
			*rowexchange_count = 0;
		}

		while (true) {
			// find pivot
			auto [ finded, finded_pivot ] = find_pivot(matrix);

			// break if not find pivot
			if ( !finded ) {
				break;
			}

			// break if submatrix only one row
			if ( matrix.rows == 1 ) {
				--matrix.rows;
				break;
			}

			// remove zero_columns
			size_t first_to_finded = std::distance(matrix.first, finded_pivot);
			size_t zerocolumn_count= first_to_finded % matrix.rowstep;
			if ( zerocolumn_count != 0 ) {
				/* remove column[0, zerocolumn_count) */
				std::advance(matrix.first, zerocolumn_count);
				matrix.cols -= zerocolumn_count;
			}

			// rowexchange first and finded, and counter
			FwdIt pivot_row = matrix.first;
			FwdIt pivot_row_last = std::next(matrix.first, matrix.cols);
			if ( pivot_row != finded_pivot ) {
				/* swap [pivot_row,pivot_row_last) and [finded_pivot, ...) */
				swap_ranges(pivot_row, pivot_row_last, finded_pivot);
				if (rowexchange_count) {
					++(*rowexchange_count);
				}
			}

			// eliminate lower rows 
			FwdIt lower_row = std::next(matrix.first, matrix.rowstep);
			for (size_t i = 1, final_i = matrix.rows-1; true; ) {
				/* [pivot_row,pivot_row_last) * -multiplier + [lower_row,...) = [lower_row, ...) */
				fma_ranges(pivot_row, pivot_row_last, -((*lower_row)/(*pivot_row)), lower_row );
				if ( i == final_i ) { 
					break; 
				}

				std::advance(lower_row, matrix.rowstep);
				++i;
			}

			// break if only column has been processed
			if ( matrix.cols == 1 ) { 
				break; 
			}
			
			// remove row[0] and remove column[0]
			std::advance(matrix.first, matrix.rowstep);
			--matrix.rows; 
			std::advance(matrix.first, 1);
			--matrix.cols;
		}

		return origin_rows - matrix.rows;
	}
	
	template<typename BidirIterator, typename Fn1, typename Fn2>
	void inplace_back_eliminate_nonsingular(MatrixArg4<BidirIterator> matrix, 
		Fn1 scale_ranges, Fn2 fma_ranges) {
		assert( matrix.cols >= matrix.rows );

		// goto back pivot
		size_t rank = matrix.rows;
		std::advance(matrix.first, (rank-1) * matrix.rowstep + (rank-1));

		while (true) {
			// test pivot
			BidirIterator pivot_row = matrix.first;
			assert( *pivot_row != 0 );
			size_t pivot_to_mid = 1 + (rank - matrix.rows);
			size_t mid_to_end   = (matrix.cols + 1) - matrix.rows - pivot_to_mid;
			BidirIterator pivot_row_mid = std::next(pivot_row, pivot_to_mid);
			BidirIterator pivot_row_end = std::next(pivot_row_mid, mid_to_end);
		
			// scale pivot row
			scale_ranges(pivot_row_mid, pivot_row_end, *pivot_row, 
				pivot_row_mid );
			*pivot_row = 1;

			// break if the only row has been scaled
			if ( matrix.rows == 1 ) {
				break;
			}
	
			// eliminate upper rows
			BidirIterator upper_row_mid = std::prev(pivot_row_mid, matrix.rowstep);
			BidirIterator upper_row     = std::prev(upper_row_mid, pivot_to_mid);
			for (size_t i = matrix.rows-1 - 1, final_i = 0; true; ) {
				fma_ranges(pivot_row_mid, pivot_row_end, -(*upper_row), 
					upper_row_mid );
				*upper_row = 0;

				if ( i == final_i ) {
					break;
				}

				upper_row_mid = std::prev(upper_row, matrix.rowstep - pivot_to_mid);
				upper_row     = std::prev(upper_row_mid, pivot_to_mid);
				--i;
			}

			// remove row[0], prev pivot
			std::advance(matrix.first, -matrix.rowstep);
			--matrix.rows;
			std::advance(matrix.first, -1);
		}
	}

	template<typename BidirIterator, typename Fn1, typename Fn2>
	void inplace_back_eliminate(MatrixArg4<BidirIterator> matrix, size_t rank, 
		Fn1 scale_ranges, Fn2 fma_ranges) {
		assert( rank <= matrix.rows );

		// ignore lower rows
		matrix.rows = rank;

		// saved origin_first and goto back pivot
		BidirIterator ori_first = matrix.first;
		std::advance(matrix.first, (rank-1) * matrix.rowstep + (rank-1));

		while ( true ) {
			// find pivot
			BidirIterator pivot_row     = matrix.first;
			BidirIterator pivot_row_end = std::next(pivot_row, (matrix.cols - matrix.rows + 1));
			size_t pivot_offset = 0;
			for ( ; pivot_row != pivot_row_end; ++pivot_row, ++pivot_offset) {
				if ( *pivot_row != 0 ) {
					break;
				}
			}

			// accelerate
			if ( pivot_offset == 0 ) {
				inplace_back_eliminate_nonsingular(
					MatrixArg4(ori_first, matrix.rows, matrix.cols, matrix.rowstep),
					scale_ranges, 
					fma_ranges
				);
				return;
			}

			// scale pivot row
			assert( *pivot_row != 0 );
			scale_ranges(pivot_row, pivot_row_end, *pivot_row, 
				pivot_row);

			// break if the only row has been scaled
			if ( matrix.rows == 1 ) {
				break;
			}

			// eliminate upper rows
			BidirIterator upper_row = std::prev(matrix.first, matrix.rowstep);
			for (size_t i = matrix.rows-1 - 1, final_i = 0; true; ) {
				fma_ranges(pivot_row, pivot_row_end, -(*upper_row), 
					upper_row);

				if ( i == final_i ) {
					break;
				}

				std::advance(upper_row, -matrix.rowstep);
				--i;
			}

			// remove row[0], prev pivot
			std::advance(matrix.first, -matrix.rowstep);
			--matrix.rows;
			std::advance(matrix.first, -1);
		}
	}

	// { eliminate source into destination }
	template<typename InputIterator, typename OutputIterator> inline
	size_t eliminate(MatrixArg4<InputIterator> source, MatrixArg2<OutputIterator> destination, size_t* rowexchange_count = nullptr) {
		bool inplace = false;
		if ( std::equality_comparable_with<InputIterator, OutputIterator> ) {
			inplace = source.first == destination.first;
		}
		
		auto _destination = MatrixArg4(
			destination.first, 
			source.rows, 
			source.cols, 
			source.rowstep
			);
		if (!inplace)
			reshape(source, _destination);
		return inplace_eliminate(_destination, rowexchange_count,
			_Find_maxabs_pivot(), _Swap_ranges(), _Fma_ranges()
			);
	}

	// { back_eliminate source into destination }
	template<typename InputIterator, typename OutputIterator> inline
	void back_eliminate(MatrixArg4<InputIterator> source, size_t rank, MatrixArg2<OutputIterator> destination) {
		bool inplace = false;
		if ( std::equality_comparable_with<InputIterator, OutputIterator> ) {
			inplace = source.first == destination.first;
		}
		
		auto _destination = MatrixArg4(
			destination.first, 
			source.rows, 
			source.cols, 
			source.rowstep
			);
		if (!inplace)
			reshape(source, _destination);
		inplace_back_eliminate(_destination, rank, 
			_Dscale_ranges(), _Fma_ranges());
	}

	// { solve source into destination }
	template<typename InputIterator, typename OutputIterator> inline
	void solve(MatrixArg4<InputIterator> source, MatrixArg2<OutputIterator> destination) {
		size_t rank = eliminate(source, destination);
		back_eliminate(destination, rank, destination);
	}
	

	template<typename FwdIt1, typename FwdIt2, typename Fn1, typename Fn2, typename Fn3>
	size_t inplace_eliminate(MatrixArg4<FwdIt1> matrix1, MatrixArg3<FwdIt2> matrix2, size_t* rowexchange_count,
		Fn1 find_pivot, Fn2 swap_ranges, Fn3 fma_ranges) {
		// saved origin_rows
		size_t origin_rows = matrix1.rows;

		// zero rowexchange count
		if (rowexchange_count) {
			*rowexchange_count = 0;
		}

		while (true) {
			// find pivot
			auto [ finded, finded_pivot] = find_pivot(matrix1);

			// break if not finded pivot
			if ( !finded ) {
				break;
			}
			
			// break if submatrix only one row
			if ( matrix1.rows == 1 ) {
				--matrix1.rows;
				break;
			}
			
			// matrix1 remove zerocolumes
			size_t first1_to_finded = std::distance(matrix1.first, finded_pivot);
			size_t zerocolumn_count = first1_to_finded % matrix1.rowstep;
			if ( zerocolumn_count != 0 ) {
				std::advance(matrix1.first, zerocolumn_count);
				matrix1.cols -= zerocolumn_count;
			}

			// rowexchange first and finded, and counter
			FwdIt1 pivot_row1     = matrix1.first;
			FwdIt1 pivot_row1_end = std::next(matrix1.first, matrix1.cols);
			FwdIt2 pivot_row2     = matrix2.first;
			FwdIt2 pivot_row2_end = std::next(matrix2.first, matrix2.cols);
			if ( pivot_row1 != finded_pivot ) {
				size_t zerorow_count = first1_to_finded / matrix1.rowstep;
				swap_ranges(pivot_row1, pivot_row1_end, finded_pivot);
				swap_ranges(pivot_row2, pivot_row2_end, std::next(pivot_row2, zerorow_count * matrix2.rowstep));
				if (rowexchange_count) {
					++(*rowexchange_count);
				}
			}

			// eliminate lower rows
			FwdIt1 lower_row1 = std::next(matrix1.first, matrix1.rowstep);
			FwdIt2 lower_row2 = std::next(matrix2.first, matrix2.rowstep);
			for (size_t i = 1, final_i = matrix1.rows-1; true; ) {
				auto neg_multiplier = -( (*lower_row1) / (*pivot_row1) );
				fma_ranges(pivot_row1, pivot_row1_end, neg_multiplier, lower_row1);
				fma_ranges(pivot_row2, pivot_row2_end, neg_multiplier, lower_row2);
				if ( i == final_i ) {
					break;
				}

				std::advance(lower_row1, matrix1.rowstep);
				std::advance(lower_row2, matrix2.rowstep);
				++i;
			}

			// break if only column has been processed
			if (matrix1.cols == 1) {
				break;
			}
			
			// remove row[0], matrix1 remove column[0]
			std::advance(matrix1.first, matrix1.rowstep);
			std::advance(matrix2.first, matrix2.rowstep);
			--matrix1.rows; 
			std::advance(matrix1.first, 1);
			--matrix1.cols;
		}

		return origin_rows - matrix1.rows;
	}

	template<typename BidIt1, typename BidIt2, typename Fn1, typename Fn2>
	void inplace_back_eliminate(MatrixArg4<BidIt1> matrix1, MatrixArg3<BidIt2> matrix2, size_t rank, 
		Fn1 scale_ranges, Fn2 fma_ranges) {
		// not singular submatrix
		assert( rank != 0);
		assert( rank <= matrix1.rows );
		assert( rank <= matrix1.cols );

		// ignore lower rows
		matrix1.rows = rank;

		// goto back pivot
		std::advance(matrix1.first, (matrix1.rows-1) * matrix1.rowstep + (rank-1));
		std::advance(matrix2.first, (matrix1.rows-1) * matrix2.rowstep);
	
		while (true) {
			// find pivot
			BidIt1 pivot_row1     = matrix1.first;
			BidIt1 pivot_row1_end = std::next(matrix1.first, matrix1.cols - matrix1.rows + 1);
			BidIt2 pivot_row2     = matrix2.first;
			BidIt2 pivot_row2_end = std::next(matrix2.first, matrix2.cols);
			for ( ; pivot_row1 != pivot_row1_end; ++pivot_row1) {
				if (*pivot_row1 != 0) {
					break;
				}
			}
			assert( *pivot_row1 != 0 );

			// scale pivot row
			const auto divisor = *pivot_row1;
			scale_ranges(pivot_row1, pivot_row1_end, divisor, pivot_row1);
			scale_ranges(pivot_row2, pivot_row2_end, divisor, pivot_row2);
			
			// break if the only row has been scaled
			if ( matrix1.rows == 1 ) {
				break;
			}
			
			// eliminate upper rows
			BidIt1 upper_row1 = std::prev(matrix1.first, matrix1.rowstep);
			BidIt2 upper_row2 = std::prev(matrix2.first, matrix2.rowstep);
			for (size_t i = matrix1.rows-1 - 1, final_i = 0; true; ) {
				const auto multiplier = *upper_row1 /*divss pivot*/;
				fma_ranges(pivot_row1, pivot_row1_end, -multiplier, upper_row1);
				fma_ranges(pivot_row2, pivot_row2_end, -multiplier, upper_row2);
				if ( i == final_i ) {
					break;
				}

				std::advance(upper_row1, -ptrdiff_t(matrix1.rowstep));
				std::advance(upper_row2, -ptrdiff_t(matrix2.rowstep));
				--i;
			}

			// remove row[0], prev pivot;
			std::advance(matrix1.first, -ptrdiff_t(matrix1.rowstep));
			std::advance(matrix2.first, -ptrdiff_t(matrix2.rowstep));
			--matrix1.rows;
			std::advance(matrix1.first, -1);
		}
	}

	// { eliminate augumented[source1|source2] into augumented[destination1|destination2] }
	template<typename InIt1, typename InIt2, typename OutIt1, typename OutIt2> inline
	size_t eliminate(MatrixArg4<InIt1> source1, MatrixArg3<InIt2> source2, MatrixArg2<OutIt1> destination1, MatrixArg2<OutIt2> destination2, size_t* rowexchange_count = nullptr) {
		bool inplace1 = false;
		bool inplace2 = false;
		if (std::equality_comparable_with<InIt1, OutIt1>) {
			inplace1 = source1.first == destination1.first;
		}
		if (std::equality_comparable_with<InIt2, OutIt2>) {
			inplace2 = source2.first == destination2.first;
		}

		auto _source2 = MatrixArg4(
			source2.first, 
			source1.rows, 
			source2.cols,
			source2.rowstep
			);
		auto _destination1 = MatrixArg4(
			destination1.first, 
			source1.rows, 
			source1.cols, 
			source1.rowstep
			);
		auto _destination2 = MatrixArg4(
			destination2.first, 
			source1.rows,
			source2.cols, 
			source2.rowstep
			);
		if (!inplace1)
			reshape(source1, _destination1);
		if (!inplace2)
			reshape(_source2, _destination2);
		return inplace_eliminate(_destination1, MatrixArg3<OutIt2>(_destination2), rowexchange_count,
			_Find_maxabs_pivot(), _Swap_ranges(), _Fma_ranges()
			);
	}

	// { back_eliminate augumented[source1|source2] into augumented[destination1|destination2] }
	template<typename InIt1, typename InIt2, typename OutIt1, typename OutIt2> inline
	void back_eliminate(MatrixArg4<InIt1> source1, MatrixArg3<InIt2> source2, size_t rank, MatrixArg2<OutIt1> destination1, MatrixArg2<OutIt2> destination2) {
		bool inplace1 = false;
		bool inplace2 = false;
		if (std::equality_comparable_with<InIt1, OutIt1>) {
			inplace1 = source1.first == destination1.first;
		}
		if (std::equality_comparable_with<InIt2, OutIt2>) {
			inplace2 = source2.first == destination2.first;
		}

		auto _source2 = MatrixArg4(
			source2.first, 
			source1.rows, 
			source2.cols,
			source2.rowstep
			);
		auto _destination1 = MatrixArg4(
			destination1.first, 
			source1.rows, 
			source1.cols, 
			source1.rowstep
			);
		auto _destination2 = MatrixArg4(
			destination2.first, 
			source1.rows, 
			source2.cols, 
			source2.rowstep
			);
		if (!inplace1)
			reshape(source1, _destination1);
		if (!inplace2)
			reshape(_source2, _destination2);
		inplace_back_eliminate(_destination1, MatrixArg3<OutIt2>(_destination2), rank,
			_Dscale_ranges(), _Fma_ranges());
	}

	// { solve augumented[source1|source2] into augumented[destination1|destination2] }
	template<typename InIt1, typename InIt2, typename OutIt1, typename OutIt2> inline
	void solve(MatrixArg4<InIt1> source1, MatrixArg3<InIt2> source2, MatrixArg2<OutIt1> destination1, MatrixArg2<OutIt2> destination2) {
		size_t rank = eliminate(source1,source2, destination1,destination2);
		back_eliminate(destination1,destination2, rank, destination1,destination2);
	}


	/*_INFORMATION({
		"Function": bareiss_solve(...),
		"Relation" : [ eliminate(...), eliminate_and_factorize(...) ],
		"Principle": {
			init pivot and divisor;
			while (size_t k = 0, last_k = diags-1, final_k = diags-2; true; ) {
				eliminate other rows using bareiss_formula;
				break if eliminated;
				next divisorand pivot;
				throw if divisor is zero;
			}
		}
	})*/

	template<typename FwdIt, typename Fn>
	void inplace_bareiss_solve(MatrixArg4<FwdIt> matrix, size_t count, Fn formula_ranges) {
		assert( matrix.rows != 0 && matrix.cols != 0 );
		assert( count < std::min(matrix.rows, matrix.cols) );
		if (count == -1) { count = std::min(matrix.rows, matrix.cols); }
		
		// init pivot and divisor
		FwdIt pivot_seek = matrix.first;
		FwdIt pivot_row = matrix.first;
		auto divisor = static_cast<typename std::iterator_traits<FwdIt>::value_type>( 1 );
		auto pivot = *pivot_seek;

		for (size_t k = 0, last_k = count-1, final_k = count-2; true; ) {
			// eliminate other rows  using bareiss_formula
			FwdIt pivot_row_last = std::next(pivot_row, matrix.cols);
			FwdIt other_row = matrix.first;
			{
				size_t final_i = matrix.rows-1;
				size_t i = 0;
				for ( ; other_row != pivot_row; ) {
					// [i, pivot)
					FwdIt multiplier_ptr = std::next(other_row, k);
					formula_ranges(other_row, pivot_row, pivot_row_last, pivot, *multiplier_ptr, divisor);
					other_row = std::next(multiplier_ptr, matrix.rowstep - k);
					++i;
				}
				// [pivot] do nothing
				while ( true ) {
					// (pivot, final)
					FwdIt multiplier_ptr = std::next(other_row, k);
					formula_ranges(other_row, pivot_row, pivot_row_last, pivot, *multiplier_ptr, divisor);
					if ( i == final_i ) {
						break;
					}
					other_row = std::next(multiplier_ptr, matrix.rowstep - k);
					++i;
				}
			}

			// break if eliminated
			if ( k == final_k ) {
				break;
			}

			// next divisor and pivot
			pivot_row = std::next(pivot_seek, matrix.rowstep - (k - 1));
			pivot_seek = std::next(pivot_row, k);
			divisor = pivot;
			pivot = *pivot_seek;
			++k;

			// throw if divisor is zero
			if ( divisor == 0 ) {
				throw std::exception("bareiss_eliminate_matrix(...)::(divisor==0)");
			}
		}
	}

	// { eliminate source into destination, count=is_determinant:diags-1|is_solve:diags }
	template<typename InputIterator, typename OutputIterator> inline
	void bareiss_solve(MatrixArg4<InputIterator> source, MatrixArg2<OutputIterator> destination, size_t count = -1) {
		bool inplace = false;
		if (std::equality_comparable_with<InputIterator, OutputIterator>) {
			inplace = source.first == destination.first;
		}

		auto _destination = MatrixArg4(
			destination.first,
			source.rows,
			source.cols,
			source.rowstep
			);
		if (!inplace)
			reshape(source, _destination);
		inplace_bareiss_solve(_destination, count,
			[]<typename FwdIt, typename Scalar>
			(FwdIt other_row, FwdIt pivot_row, FwdIt pivot_row_end, const Scalar& pivot, const Scalar& multiplier, const Scalar& divisor) {
				std::transform(pivot_row, pivot_row_end, other_row,
					other_row,
					[pivot, multiplier, divisor](const Scalar& pivot_row_element, const Scalar& other_row_element) {
						return (other_row_element * pivot - pivot_row_element * multiplier) / divisor;
					}
				);
			}
			);
	}

	// parallel bareiss elimination ...



	/*_INFORMATION({
		"Function": eliminate_and_factorize(...),
		"Relation": [ eliminate(...), bareiss_eliminate(...) ],
		"Principle": {
			saved origin_rows and setup inplace;
			(Swapable)init rowexchange and zero rowexchange_count;
			while (true) {
				find pivot; (Swapable)
				break if not finded pivot; 
				break if has zerocolumn; (Swapable)
				break if submatrix only one row;
				rowexchange first and finded, and counter; (Swapable)
				eliminate lower rows and store multipliers;
				break if only column has been processed;
				remove row[0], remove column[0]
			}
		}
	})*/

	template<typename FwdIt1, typename FwdIt2, typename Fn>
	size_t inplace_eliminate_and_factorize(MatrixArg4<FwdIt1> matrix, MatrixArg2<FwdIt2> multipliers, 
		Fn fma_ranges) {
		// saved origin_rows and setup inplace
		size_t origin_rows = matrix.rows;
		bool inplace = false;
		if constexpr (std::equality_comparable_with<FwdIt1, FwdIt2>) {
			inplace = (matrix.first == multipliers.first);
		}
		
		while ( true ) {
			// break if error pivot
			if ( *matrix.first == 0 ) {
				break;
			}

			// break if submatrix only one row
			if ( matrix.rows == 1 ) {
				if (!inplace) { *multipliers.first = 1; }
				--matrix.rows;
				break;
			}
			
			// eliminate lower rows and store multipliers
			if ( !inplace ) { *multipliers.first = 1; }
			FwdIt1 pivot_row        = matrix.first;
			FwdIt1 pivot_row_last   = std::next(matrix.first, matrix.cols);
			FwdIt1 lower_row        = std::next(matrix.first, matrix.rowstep);
			FwdIt2 lower_multiplier = std::next(multipliers.first, multipliers.rowstep);
			for(size_t i = 1, final_i = matrix.rows-1; true; ) { 
				auto multiplier = *lower_row / *pivot_row;// !!!inplace requires
				fma_ranges(pivot_row, pivot_row_last, -multiplier, lower_row);
				*lower_multiplier = multiplier;
				if ( i == final_i ) {
					break;
				}

				std::advance(lower_row, matrix.rowstep);
				std::advance(lower_multiplier, multipliers.rowstep);
				++i;
			}

			// break if only column has been processed
			if ( matrix.cols == 1 ) {
				break;
			}

			// remove row[0], remove column[0]
			std::advance(matrix.first, matrix.rowstep);
			std::advance(multipliers.first, multipliers.rowstep);
			--matrix.rows;
			std::advance(matrix.first, 1);
			std::advance(multipliers.first, 1);
			--matrix.cols;
		}

		return origin_rows - matrix.rows;
	}
	
	// { factorize: multipliers * matrix = origin_matrix }
	template<typename InputIterator, typename OutputIterator, typename OutIt2> inline
	size_t eliminate_and_factorize(MatrixArg4<InputIterator> source, MatrixArg2<OutputIterator> destination, MatrixArg2<OutIt2> multipliers) {
		bool inplace = false;
		if (std::equality_comparable_with<InputIterator, OutputIterator>) {
			inplace = source.first == destination.first;
		}

		auto _destination = MatrixArg4(
			destination.first,
			source.rows,
			source.cols,
			source.rowstep
			);
		if (!inplace) {
			reshape(source, _destination);
			}
		return inplace_eliminate_and_factorize(_destination, multipliers,
			_Fma_ranges()
			);
	}


	template<typename FwdIt1, typename FwdIt2, typename FwdIt3, typename Fn1, typename Fn2, typename Fn3>
	size_t inplace_eliminate_and_factorize(MatrixArg4<FwdIt1> matrix, MatrixArg2<FwdIt2> multipliers, MatrixArg2<FwdIt3> rowexchange, size_t* rowexchange_count, 
		Fn1 find_pivot, Fn2 swap_ranges, Fn3 fma_ranges) {
		// saved ori_diags and setup inplace
		size_t origin_rows = matrix.rows;
		bool inplace = false;
		if constexpr (std::equality_comparable_with<FwdIt1, FwdIt2>) {
			inplace = (matrix.first == multipliers.first);
		}
		
		// init rowexchange and zero rowexchange_count
		set_identity(MatrixArg4(rowexchange.first, matrix.rows, matrix.cols, rowexchange.rowstep));
		if (rowexchange_count) {
			*rowexchange_count = 0;
		}
		
		FwdIt2 multipliers_first_row = multipliers.first;
		while ( true ) {
			// find pivot
			auto [ finded, finded_pivot ] = find_pivot(matrix);

			// break if not finded pivot
			if ( !finded ) {
				break;
			}
			
			// break if has zerocolumn
			size_t first_to_finded = std::distance(matrix.first, finded_pivot);
			if ( (first_to_finded % matrix.rowstep) != 0 ) {
				break;
			}

			// break if submatrix only one row
			if ( matrix.rows == 1 ) {
				if ( !inplace ) { *multipliers.first = 1; }
				--matrix.rows;
				break;
			}

			// rowexchange first and finded, and counter
			FwdIt1 pivot_row      = matrix.first;
			FwdIt1 pivot_row_last = std::next(matrix.first, matrix.cols);
			if ( pivot_row != finded_pivot ) {
				size_t zerorow_count = first_to_finded / matrix.rowstep;
				swap_ranges(pivot_row, pivot_row_last, finded_pivot);
				swap_ranges(multipliers_first_row, multipliers.first, 
							std::next(multipliers_first_row, zerorow_count * multipliers.rowstep));
				{// swap_ranges(permute_matrix)
					FwdIt3 leftr_first = std::next(rowexchange.first, 0 * rowexchange.rowstep);
					FwdIt3 rightr_first = std::next(rowexchange.first, zerorow_count * rowexchange.rowstep);

					FwdIt3 leftr_source = leftr_first;
					FwdIt3 leftr_destination = rightr_first;
					while (true) {
						if ( *leftr_source != 0 ) {
							break;
						}

						std::advance(leftr_source, 1);
						std::advance(leftr_destination, 1);
					}

					FwdIt3 rightr_source = rightr_first;
					FwdIt3 rightr_destination = leftr_first;
					while (true) {
						if (*rightr_source != 0) {
							break;
						}

						rightr_source = std::next(rightr_source);
						rightr_destination = std::next(rightr_destination);
					}

					std::swap(*leftr_source, *leftr_destination);
					std::swap(*rightr_source, *rightr_destination);
				}
				if (rowexchange_count) {
					++(*rowexchange_count);
				}
			}

			// eliminate lower rows and store multipliers
			if ( !inplace ) { *multipliers.first = 1; }
			FwdIt1 lower_row        = std::next(matrix.first, matrix.rowstep);
			FwdIt2 lower_multiplier = std::next(multipliers.first, multipliers.rowstep);
			for(size_t i = 1, final_i = matrix.rows-1; true; ) { 
				auto multiplier = *lower_row / *pivot_row;// !!!inplace requires
				fma_ranges(pivot_row, pivot_row_last, -multiplier, lower_row);
				*lower_multiplier = multiplier;
				if ( i == final_i ) {
					break;
				}

				std::advance(lower_row, matrix.rowstep);
				std::advance(lower_multiplier, multipliers.rowstep);
				++i;
			}

			// break if only column has been processed;
			if ( matrix.cols == 1 ) {
				break;
			}

			// remove row[0], remove column[0]
			std::advance(matrix.first, matrix.rowstep);
			std::advance(rowexchange.first, rowexchange.rowstep);
			std::advance(multipliers_first_row, multipliers.rowstep);
			--matrix.rows;
			std::advance(matrix.first, 1);
			multipliers.first = std::next(multipliers_first_row, origin_rows - matrix.rows);
			--matrix.cols;
		}

		return origin_rows - matrix.rows;
	}
	
	// { factorize: multipliers * matrix = rowexchange * origin_matrix }
	template<typename InputIterator, typename OutputIterator, typename OutIt2, typename OutIt3> inline
	size_t eliminate_and_factorize(MatrixArg4<InputIterator> source, MatrixArg2<OutputIterator> destination, MatrixArg2<OutIt2> multipliers, MatrixArg2<OutIt3> rowexchange, size_t* rowexchange_count = nullptr) {
		bool inplace = false;
		if (std::equality_comparable_with<InputIterator, OutputIterator>) {
			inplace = source.first == destination.first;
		}

		auto _destination = MatrixArg4(
			destination.first,
			source.rows,
			source.cols,
			source.rowstep
			);
		if (!inplace) {
			reshape(source, _destination);
			}
		return inplace_eliminate_and_factorize(_destination, 
			multipliers, rowexchange, rowexchange_count,
			_Find_maxabs_pivot(), _Swap_ranges(), _Fma_ranges()
			);
	}



	/*_INFORMATION({
		"Function": determinant(...),
		"Relation": [ bareiss_eliminate(...), eliminate(...), "Cramer-theorem" ],
		"Principle": {
			det = sum( pow(-1,count_invseq) * product(matrix.at(i,seq[i])) );
		}
	})*/

	template<typename FwdIt>
	auto make_permutations(FwdIt first, FwdIt last)
		-> std::vector<std::vector<typename std::iterator_traits<FwdIt>::value_type>>
	{
		using T = typename std::iterator_traits<FwdIt>::value_type;
		size_t n = std::distance(first, last);

		std::vector<std::vector<T>> sequences;
		sequences.resize(fact<size_t>(n));

		// copy base_sequence
		for ( ; first != last; ++first) {
			sequences[0].push_back(*first);
		}

		size_t subsequences = 1;
		size_t subsequence_size = n;
		size_t subsequence_stride = sequences.size() / subsequence_size;
		while (true) {
			// foreach subsequence
			for (size_t i = 0; i != subsequences; ++i) {
				size_t subsequence_offset = i * subsequence_size * subsequence_stride;

				auto& subsequence_base = sequences[subsequence_offset];

				for (size_t k = 1; k != subsequence_size; ++k) {
					// subsequence[i] = subsequence_base swaped from [n-subsequence_size, ...) 
					auto& the_sequence = sequences[subsequence_offset + k * subsequence_stride];
					the_sequence = subsequence_base;
					for (size_t j = k; j != 0; --j) {
						std::swap(the_sequence[(n-subsequence_size) + j],
									the_sequence[(n-subsequence_size) + j-1]);
					}
				}
			}

			// next level subsequence
			subsequence_size -= 1;
			// break if subsequence only one
			if (subsequence_size == 1) {
				break;
			}
			subsequence_stride = subsequence_stride / subsequence_size; /*fact<size_t>(subsequence_size - 1);*/
			subsequences = sequences.size() / (subsequence_size * subsequence_stride);
		}

		return std::move(sequences);
	}

	template<typename FwdIt>
	auto determinant(FwdIt first, size_t diags, size_t rowstep)
			-> typename std::iterator_traits<FwdIt>::value_type
		{
			// make base sequence
			auto base_sequence = std::vector<size_t>(diags);
			for (size_t i = 0; i != diags; ++i) {
				base_sequence[i] = i;
			}
		
			// make full sequences
			auto sequences = make_permutations(base_sequence.begin(), base_sequence.end());
		
			// det = pow(-1,count_invseq) * product(matrix[i++][seq[i++]])
			auto det = static_cast<typename std::iterator_traits<FwdIt>::value_type>( 0 );
			for (size_t k = 0; k != sequences.size(); ++k) {
				auto product = static_cast<decltype(det)>(1);
				auto       seq      = sequences[k].begin();
				const auto seq_end  = sequences[k].end();
				FwdIt      each_row = first;
				while (true) {
					product *= *std::next(each_row, *seq);
					if (++seq == seq_end) {
						break;
					}
					each_row = std::next(each_row, rowstep);
				}

				uint32_t count_invseq = 0;
				for (auto it = std::is_sorted_until(sequences[k].begin(), sequences[k].end()); 
					it != sequences[k].end(); ++it) 
				{
					count_invseq += std::count_if(sequences[k].begin(), it, 
						[it](size_t s) { return (s > (*it)); });
				}

				det += product * pow(-1, count_invseq);
			}

			return det;
		}



#ifdef _VECTOR_
	template<typename InputIterator1, typename InputIterator2, typename _InIt3, typename OutputIterator,
		typename _Fn1, typename _Fn2>
	void diagmatrix_fixed_point_iterate(
		InputIterator1 F_begin, size_t F_diags, size_t F_rowstep, 
		OutputIterator X_begin, size_t X_cols, size_t X_rowstep,
		InputIterator2 Y_begin, size_t Y_rowstep,
		_InIt3 X0_begin, size_t X0_rowstep,
		size_t step,
		_Fn1 _r_assign_mulab,
		_Fn2 _r_addassign_mulab) {

		using _Ty = decltype(*F_begin + *F_begin);

		// copy X, X0
		OutputIterator X_ptr  = X_begin;
		_InIt3 X0_ptr = X0_begin;
		for (size_t i = 0; true; ++i) {
			std::copy(X0_ptr, std::next(X0_ptr, X_cols), X_ptr);

			if (i == (F_diags - 1)) {
				break;
			}
			X_ptr = std::next(X_ptr, X_rowstep);
			X0_ptr = std::next(X0_ptr, X0_rowstep);
		}

		// register
		auto R = std::vector<double>(F_diags * X_cols);

		for (size_t k = 0; k != step; ++k) {

			// R = invD * (Y - (L+U)*X);
			InputIterator1 F_row = F_begin;
			OutputIterator X_row = X_begin;
			InputIterator2 Y_row = Y_begin;
			_Ty*   R_row = R.data();
			for (size_t i = 0; true; ++i) {
				// 1. decompose LDU
				InputIterator1 L_begin = F_row;
				InputIterator1 L_end   = std::next(F_row, i);
				auto   D       = *L_end;
				InputIterator1 U_begin = std::next(L_end);
				InputIterator1 U_end   = std::next(L_end, F_diags - i);
				// 2. setup
				_Ty* R_row_d = R_row + X_cols;

				// 3. (L+U)*X
				OutputIterator X_ptr = X_begin;
				if (L_begin != L_end) {
					// L*x
					_r_assign_mulab( R_row, R_row_d, X_ptr, *L_begin );
					L_begin = std::next(L_begin);
					X_ptr   = std::next(X_ptr, X_rowstep);
					while (L_begin != L_end) {
						_r_addassign_mulab( R_row, R_row_d, X_ptr, *L_begin );
						L_begin = std::next(L_begin);
						X_ptr   = std::next(X_ptr, X_rowstep);// safe
					}
					// D*x
					X_ptr = std::next(X_ptr, U_begin != U_end ? X_rowstep : 0);
					// U*x...
				} else {
					// L*x...
					// D*x
					X_ptr = std::next(X_ptr, X_rowstep);
					// U*x
					_r_assign_mulab( R_row, R_row_d, X_ptr, *U_begin );
					U_begin = std::next(U_begin);
					X_ptr   = std::next(X_ptr, U_begin != U_end ? X_rowstep : 0);
				}
				// -------------------
				while (U_begin != U_end) {
					_r_addassign_mulab( R_row, R_row_d, X_ptr, *U_begin );
					U_begin = std::next(U_begin);
					X_ptr   = std::next(X_ptr, U_begin != U_end ? X_rowstep : 0);
				}

				// 4. Y - (L+U)*X
				std::transform( R_row, R_row_d, Y_row, R_row,
					[](_Ty lux, _Ty y) { return y - lux; } );

				// 5. invD * (Y - (L+U)*X)
				std::transform( R_row, R_row_d, R_row,
					[D](_Ty y_lux) { return y_lux / D; } );

				// 6. next
				if (i == (F_diags - 1)) {
					break;
				}
				F_row = std::next(F_row, F_rowstep);
				X_row = std::next(X_row, X_rowstep);
				Y_row = std::next(Y_row, Y_rowstep);
				R_row = std::next(R_row, X_cols);
			}

			// X = R
			OutputIterator X_ptr = X_begin;
			_InIt3 R_ptr = R.data();
			for (size_t i = 0; true; ++i) {
				std::copy(R_ptr, std::next(R_ptr, X_cols), X_ptr);

				if (i == (F_diags - 1)) {
					break;
				}
				X_ptr = std::next(X_ptr, X_rowstep);
				R_ptr = std::next(R_ptr, X_cols);
			}

		}
	}

	template<typename InputIterator1, typename InputIterator2, typename _InIt3, typename OutputIterator>
	void diagmatrix_fixed_point_iterate(
		InputIterator1 F_begin, size_t F_diags, size_t F_rowstep,
		OutputIterator X_begin, size_t X_cols, size_t X_rowstep,
		InputIterator2 Y_begin, size_t Y_rowstep,
		_InIt3 X0_begin, size_t X0_rowstep,
		size_t step) {

		auto _r_assign_mulab = 
			[](auto r_begin, auto r_end, InputIterator1 b_begin, auto&& a) {
				std::transform(
					r_begin, r_end, 
					b_begin, 
					r_begin,
					[a](auto&& r, auto&& b) { return a * b; });
			};

		auto _r_addassign_mulab =
			[](auto r_begin, auto r_end, InputIterator1 b_begin, auto&& a) {
				std::transform(
					r_begin, r_end, 
					b_begin, 
					r_begin,
					[a](auto&& r, auto&& b) { return a * b + r; });
			};

		diagmatrix_fixed_point_iterate(
			F_begin, F_diags, F_rowstep,
			X_begin, X_cols, X_rowstep,
			Y_begin, Y_rowstep,
			X0_begin, X0_rowstep,
			step,
			_r_assign_mulab, 
			_r_addassign_mulab);
	}
#endif// #include <vector>

	// { Cholesky factorization }

	// { QR factorization, Least squares }

	// { Conjugate Gradient iterate }








	// undeterminant
	
	template<matrixconcept genMatrix, typename Iter>
	genMatrix make_vandermonde_matrix(Iter x_itr, size_t x_size) {
		using     genReal = matrix_scalar_t<genMatrix>;
		genMatrix matrix  = genMatrix();
		size_t    rowstep = matrix.rowstep();
		auto      itr     = matrix.begin();
		for (size_t i = 0; true; ) {
			genReal x = static_cast<genReal>(*x_itr++);
			for (size_t j = 0; j != x_size; ++j) {
				*itr++ = pow(x, static_cast<genReal>(j));
			}

			if (++i == x_size) {
				break;
			}
			
			itr = std::next(itr, rowstep - x_size);
		}

		return std::move(matrix);
	}

	template<typename OutputIterator, typename InIt1, typename InIt2>
	void set_remap_matrix(
		OutputIterator scalings, size_t rowstep,
		InIt1 lowers, 
		InIt1 uppers, 
		InIt2 newlowers, 
		InIt2 newuppers,
		size_t equations
	) {
		OutputIterator translations = std::next(scalings, equations);
		while (true) {
			// (x-lower)/(upper-lower)*(newupper-newlower) + newlower
			// (x-lower)*ratio + newlower
			// x*ratio - lower*ratio + newlower
			auto lower    = *lowers;
			auto newlower = *newlowers;
			auto ratio    = (*newuppers - *newlowers)/(*uppers - *lowers);
			*scalings     = ratio;
			*translations = -lower*ratio + newlower;

			// break if equations is empty
			if (--equations == 0) {
				break;
			}

			// next equation next variable
			scalings = std::next(scalings, rowstep + 1);
			++lowers; ++uppers;
			++newlowers; ++newuppers;
			// next equation
			translations = std::next(translations, rowstep);
		}

		// invariable
		*std::next(translations,rowstep) = 1;
	}

	template<matrixconcept genMatrix, typename InIt1, typename InIt2>
	void set_remap_matrix(genMatrix& matrix, InIt1 lowers, InIt1 uppers, InIt2 newlowers, InIt2 newuppers, size_t equations) {
		assert(matrix.rows() >= equations+1);
		assert(matrix.cols() >= equations+1);
		set_remap_matrix(matrix.begin(), matrix.rowstep(),
			lowers, uppers, 
			newlowers, newuppers, 
			equations
		);
	}

#ifndef make_matrix3x3
#define make_matrix3x3(MType, S, m00, m01, m02, m10, m11, m12, m20, m21, m22) \
	MType{ static_cast<S>(m00), static_cast<S>(m01), static_cast<S>(m02), \
		   static_cast<S>(m10), static_cast<S>(m11), static_cast<S>(m12), \
		   static_cast<S>(m20), static_cast<S>(m21), static_cast<S>(m22) } 
#endif

#ifndef make_matrix4x4
#define make_matrix4x4(MType, S, m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33) \
	MType{ static_cast<S>(m00), static_cast<S>(m01), static_cast<S>(m02), static_cast<S>(m03),  \
		   static_cast<S>(m10), static_cast<S>(m11), static_cast<S>(m12), static_cast<S>(m13),  \
		   static_cast<S>(m20), static_cast<S>(m21), static_cast<S>(m22), static_cast<S>(m23),  \
		   static_cast<S>(m30), static_cast<S>(m31), static_cast<S>(m32), static_cast<S>(m33) } 
#endif

}// namespace calculation
/*
// { diagmatrixA{...} * matrixB(DxP) = matrixR(DxP), inplace }
template<typename FwdIt1, typename FwdIt2, typename FwdIt3, typename Fn>
void linear_combination(
	FwdIt1 diagmatrixA_itr, 
	size_t diagmatrixA_diags,
	FwdIt2 matrixB_itr, 
	size_t matrixB_cols, 
	size_t matrixB_rowstep,
	FwdIt3 matrixR_itr, 
	size_t matrixR_rowstep, Fn transform_operator
) {
	FwdIt1 diagmatrixA_back = std::next(diagmatrixA_itr, diagmatrixA_diags - 1);
	while (true) {
		// Rvector = weight * vector
		transform_operator(*diagmatrixA_itr, matrixB_itr, 
			matrixR_itr, std::next(matrixR_itr, matrixB_cols) );
		// cmp
		break_if_true(diagmatrixA_itr == diagmatrixA_back);
		// next major
		diagmatrixA_itr = std::next(diagmatrixA_itr);
		// next row
		matrixB_itr = std::next(matrixB_itr, matrixB_rowstep);
		matrixR_itr = std::next(matrixB_itr, matrixR_rowstep);
		// jump
	}
}

// { diagmatrixA{...} * matrixB(DxP) = matrixR(DxP), inplace }
template<typename FwdIt1, typename FwdIt2, typename FwdIt3>
void linear_combination(
	FwdIt1 diagmatrixA_itr, 
	size_t diagmatrixA_diags,
	FwdIt2 matrixB_itr, 
	size_t matrixB_cols, 
	size_t matrixB_rowstep, 
	FwdIt3 matrixR_itr, 
	size_t matrixR_rowstep
) {
	auto transform_operator = [](FwdIt2 vector_itr, auto&& weight, FwdIt3 result_itr, FwdIt3 result_itr_end) {
		std::transform(
			result_itr, result_itr_end, 
			vector_itr,
			result_itr, 
			[weight](auto&& r, auto&& b) { return weight * b; } 
		);
	};
	linear_combination(
		diagmatrixA_itr, 
		diagmatrixA_diags, 
		matrixB_itr, 
		matrixB_cols, 
		matrixB_rowstep, 
		matrixR_itr, 
		matrixR_rowstep, 
		transform_operator
	);
}
*/