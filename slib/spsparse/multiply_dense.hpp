
/** @defgroup multiply_dense multiply_dense.hpp
@brief Multiplication of sparse matrices by dense vectors.

@{
*/

#if 0

/** (Sparse Matrix) * (Dense Vector) */
template<class VectorCooMatrixT, class AccumulatorT>
void multiply(
	VectorCooMatrixT const &M,
	blitz::Array<double,1> const &x,
	AccumulatorT &y,
	bool handle_nan = false,
	bool transpose = false)
{
	for (auto ii = M.begin(); ii != M.end(); ++ii) {
		std::array<IndexT,2> index(ii.index());
		if (transpose) std::swap(index[0], index[1]);
		double val = ii.val() * x(index[1]);
		if (!handle_nan || !(std::isnan(val) || std::isinf(val))) {
			y.add(  {{index[0]}},  val);
		}
	}
}

template<class VectorCooMatrixT, class AccumulatorT>
void multiplyT(
	VectorCooMatrixT const &M,
	blitz::Array<double,1> const &x,
	AccumulatorT &y,
	bool handle_nan = false)
{ return multiply(M, x, y, handle_nan, true); }
#endif

// ---------------------------------------------------

/** @} */
