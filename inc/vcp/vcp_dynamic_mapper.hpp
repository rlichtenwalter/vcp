/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_VCP_DYNAMIC_MAPPER
#define VCP_VCP_DYNAMIC_MAPPER

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include <vcp/multirelational_graph.hpp>
#include <vcp/square_matrix.hpp>

namespace vcp {
	
template <std::size_t base,std::size_t exponent>
struct TMP_power {
	enum {
		value_matrix = base * TMP_power<base,exponent - 1>::value_matrix
	};
};

template <std::size_t base>
struct TMP_power<base,0> {
	enum {
		value_matrix = 1
	};
};

template <std::size_t value1, std::size_t value2>
struct TMP_min {
	enum {
		value_matrix = value1 < value2 ? value1 : value2
	};
};

template <std::size_t value1, std::size_t value2>
struct TMP_max {
	enum {
		value_matrix = value1 > value2 ? value1 : value2
	};
};

template <std::size_t n,std::size_t r,bool d>
class vcp_dynamic_mapper {
	public:
		typedef typename multirelational_graph<r>::connectivity_address_type connectivity_address_type;
		typedef typename std::conditional<
				n*(n-1)*r*(d+1)/2 <= CHAR_BIT*sizeof(std::size_t),
				std::size_t,
				boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
						n*(n-1)*r*(d+1)/2,
						n*(n-1)*r*(d+1)/2,
						boost::multiprecision::unsigned_magnitude,
						boost::multiprecision::unchecked,
						void> >
				>::type subgraph_address_type;
		vcp_dynamic_mapper();
		constexpr subgraph_address_type subgraph_count() const;
		subgraph_address_type subgraph_address( square_matrix<connectivity_address_type,n> const & connectivity ) const;
		subgraph_address_type canonical_subgraph_address( square_matrix<connectivity_address_type,n> const & connectivity ) const;
		square_matrix<connectivity_address_type,n> element_structure( subgraph_address_type const & address ) const;
	private:
		square_matrix<std::size_t,n> value_matrix;
		subgraph_address_type canonical_subgraph_address( square_matrix<connectivity_address_type,n> const & connectivity, subgraph_address_type subgraph_address ) const;
};

template <std::size_t n,std::size_t r,bool d>
vcp_dynamic_mapper<n,r,d>::vcp_dynamic_mapper() {
	std::size_t index( 0 );
	for( std::size_t row( 0 ); row < n; ++row ) {
		value_matrix( row, row ) = 0;
		for( std::size_t column( row + 1 ); column < n; ++column ) {
			value_matrix( row, column ) = r * index++;
			value_matrix( column, row ) = d ? (r * index++) : value_matrix( row, column );
		}
	}
}

template <std::size_t n,std::size_t r,bool d>
constexpr typename vcp_dynamic_mapper<n,r,d>::subgraph_address_type vcp_dynamic_mapper<n,r,d>::subgraph_count() const {
	return subgraph_address_type(1) << n*(n-1)*r*(d+1)/2;
}


template <std::size_t n,std::size_t r,bool d>
typename vcp_dynamic_mapper<n,r,d>::subgraph_address_type vcp_dynamic_mapper<n,r,d>::subgraph_address( square_matrix<connectivity_address_type,n> const & connectivity ) const {
	subgraph_address_type subgraph_address( 0 );
	for( std::size_t row( 0 ); row < n; ++row ) {
		for( std::size_t column( d ? 0 : row + 1 ); column < n; ++column ) {
			if( row != column ) {
				subgraph_address += subgraph_address_type( connectivity( row, column ) ) << value_matrix( row, column );
			}
		}
	}
	return subgraph_address;
}

template <std::size_t n,std::size_t r,bool d>
typename vcp_dynamic_mapper<n,r,d>::subgraph_address_type vcp_dynamic_mapper<n,r,d>::canonical_subgraph_address( square_matrix<connectivity_address_type,n> const & connectivity ) const {
	return canonical_subgraph_address( connectivity, subgraph_address( connectivity ) );
}

template <std::size_t n,std::size_t r,bool d>
typename vcp_dynamic_mapper<n,r,d>::subgraph_address_type vcp_dynamic_mapper<n,r,d>::canonical_subgraph_address( square_matrix<connectivity_address_type,n> const & connectivity, subgraph_address_type subgraph_address ) const {
	std::array<std::size_t,n> permuter;
	for( std::size_t row( 0 ); row < n; ++row ) {
		permuter[row] = row;
	}
	while( std::next_permutation( permuter.begin() + 2, permuter.end() ) ) {
		subgraph_address_type isomorphism_address( 0 );
		for( std::size_t row( 0 ); row < n; ++row ) {
			for( std::size_t column( 0 ); column < n; ++column ) {
				if( row != column ) {
					isomorphism_address += subgraph_address_type( connectivity( row, column ) ) << value_matrix( permuter[row], permuter[column] );
				}
			}
		}
		subgraph_address = std::min( subgraph_address, isomorphism_address );
	}
	return subgraph_address;
}

template <std::size_t n,std::size_t r,bool d>
square_matrix<typename vcp_dynamic_mapper<n,r,d>::connectivity_address_type,n> vcp_dynamic_mapper<n,r,d>::element_structure( subgraph_address_type const & address ) const {
	subgraph_address_type r_pset( vcp_dynamic_mapper::subgraph_address_type(1) << r );
	square_matrix<connectivity_address_type,n> matrix;
	std::size_t row( 0 );
	std::size_t column( 1 );
	while( address > 0 ) {
		matrix( row, column ) = address & (r_pset - 1);
		address >>= r_pset;
		if( d ) {
			std::swap( row, column );
			if( row < column ) {
				++column;
				if( column >= n ) {
					++row;
					column = row + 1;
				}
			}
		} else {
			++column;
			if( column >= n ) {
				++row;
				column = row + 1;
			}
		}
	}
	return matrix;
}

}

#endif

