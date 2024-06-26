// Copyright 2020-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_LAPACK_SYEV_HPP
#define BOOST_MULTI_ADAPTORS_LAPACK_SYEV_HPP

#include <boost/multi/adaptors/blas/filling.hpp"
#include <boost/multi/adaptors/lapack/core.hpp"

#include <boost/multi/config/NODISCARD.hpp"

#include <cassert>

namespace boost {
namespace multi {
namespace lapack {

using blas::filling;

using ::core::syev;

template<class Array2D, class Array1D, class Array1DW>
auto syev(blas::filling uplo, Array2D&& a, Array1D&& w, Array1DW&& work)
	-> decltype(syev('V', uplo == blas::filling::upper ? 'L' : 'U', size(a), base(a), stride(a), base(w), base(work), size(work), std::declval<int&>()), a({0L, 1L}, {0L, 1L})) {
	assert(size(work) >= std::max(1L, 3 * size(a) - 1L));
	assert(size(a) == size(w));
	assert(stride(w) == 1);
	assert(stride(work) == 1);

	if(size(a) == 0)
		return std::forward<Array2D>(a)();

	int info = -1;

	if(stride(rotated(a)) == 1) {
		syev('V', uplo == blas::filling::upper ? 'L' : 'U', size(a), base(a), stride(a), base(w), base(work), size(work), info);
	} else if(stride(a) == 1) {
		syev('V', uplo == blas::filling::upper ? 'U' : 'L', size(a), base(a), stride(rotated(a)), base(w), base(work), size(work), info);
	} else {
		assert(0);
	}  // case not contemplated by lapack

	if(info < 0) {
		assert(0);
	}  // bad argument

	return std::forward<Array2D>(a)({0, size(a) - info}, {0, size(a) - info});
}

template<class Array2D, class Array1D, class Array1DW = typename std::decay_t<Array1D>::decay_type>
auto syev(blas::filling uplo, Array2D&& a, Array1D&& w)
	-> decltype(syev(uplo, std::forward<Array2D>(a), std::forward<Array1D>(w), Array1DW(std::max(1L, 3 * size(a) - 1L), get_allocator(w)))) {
	return syev(uplo, std::forward<Array2D>(a), std::forward<Array1D>(w), Array1DW(std::max(1L, 3 * size(a) - 1L), get_allocator(w)));
}  // TODO(correaa) obtain automatic size from lapack info routine

template<class Array2D, class Array1D>
NODISCARD("because input array is const, output gives eigenvectors")
typename Array2D::decay_type syev(blas::filling uplo, Array2D const& a, Array1D&& w) {
	auto ret = a.decay();
	auto l   = syev(uplo, ret, std::forward<Array1D>(w));
	if(size(l) != size(a))
		assert(0);  // failed
	return ret;
}

template<class Array2D>
NODISCARD("because input array is const, output gives eigenvalues")
auto syev(blas::filling uplo, Array2D&& a) {
	multi::array<typename std::decay_t<Array2D>::element_type, 1, decltype(get_allocator(a))> eigenvalues(size(a), get_allocator(a));
	syev(uplo, std::forward<Array2D>(a), eigenvalues);
	return eigenvalues;
}

template<class Array2D>
NODISCARD("because input array is const, output gives a structured binding of eigenvectors and eigenvactor")
auto syev(blas::filling uplo, Array2D const& a) {
	struct {
		typename Array2D::decay_type eigenvectors;
		typename Array2D::value_type eigenvalues;
	} ret{a, typename Array2D::value_type(size(a), get_allocator(a))};
	auto&& l = syev(uplo, ret.eigenvectors, ret.eigenvalues);
	assert(size(l) == size(a));
	return ret;
}

}  // namespace lapack
}  // namespace multi
}  // namespace boost
#endif  // BOOST_MULTI_ADAPTORS_LAPACK_SYEV_HPP
