/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_VCP_STATIC_MAPPER
#define VCP_VCP_STATIC_MAPPER

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <vector>
#include <vcp/square_matrix.hpp>

namespace vcp {

class vcp_static_mapper {
	public:
		static std::size_t subgraph_count( std::size_t n, std::size_t r, bool d );
		vcp_static_mapper( std::size_t n, std::size_t r, bool d );
		std::size_t n() const;
		std::size_t r() const;
		bool d() const;
		template <std::size_t n> std::size_t subgraph_address( square_matrix<std::size_t,n> const & connectivity ) const;
		template <std::size_t n> std::size_t element_address( square_matrix<std::size_t,n> const & connectivity ) const;
		std::size_t subgraph_address( square_matrix<std::size_t> const & connectivity ) const;
		std::size_t element_address( square_matrix<std::size_t> const & connectivity ) const;
		std::size_t element_address( std::size_t subgraph_address ) const;
		square_matrix<std::size_t> element_structure( std::size_t element_address ) const;
		friend std::ostream & operator<<( std::ostream & os, vcp_static_mapper const & mapper );
	private:
		std::size_t n_;
		std::size_t r_;
		bool d_;
		std::vector<std::size_t> map;
		square_matrix<std::size_t> value_matrix;
};

std::size_t vcp_static_mapper::subgraph_count( std::size_t n, std::size_t r, bool d ) {
	return std::pow( 2, n*(n-1)*r*(d+1)/2 );
}

vcp_static_mapper::vcp_static_mapper( std::size_t n, std::size_t r, bool d ) : n_(n), r_(r), d_(d), map( vcp_static_mapper::subgraph_count(n,r,d) ), value_matrix( n_ )  {
	std::size_t r_pset( std::pow( 2, r_ ) );
	std::size_t index( 0 );
	for( std::size_t row( 0 ); row < n_; ++row ) {
		value_matrix( row, row ) = 0;
		for( std::size_t column( row + 1 ); column < n_; ++column ) {
			value_matrix( row, column ) = std::pow( r_pset, index++ );
			value_matrix( column, row ) = d_ ? std::pow( r_pset, index++ ) : value_matrix( row, column );
		}
	}

	square_matrix<std::size_t> connectivity( n_ );
	std::vector<std::size_t> permuter( n_ );	
	
	std::size_t element_address( 0 );
	while( !connectivity( 0, 0 ) ) {
		std::size_t subgraph_address( 0 );
		for( std::size_t row( 0 ); row < n_; ++row ) {
			for( std::size_t column( d_ ? 0 : row + 1 ); column < n_; ++column ) {
				subgraph_address += connectivity( row, column ) * value_matrix( row, column );
			}
		}
		if( map[ subgraph_address ] == 0 ) { // we have not encountered this subgraph or its isomorphic equivalents before
			map[ subgraph_address ] = element_address;
			for( std::size_t i( 0 ); i < n_; ++i ) { // reset the permuter vector
				permuter[i] = i;
			}
			while( std::next_permutation( permuter.begin() + 2, permuter.end() ) ) { // find all isomorphisms
				std::size_t isomorphism_address( 0 );
				for( std::size_t row( 0 ); row < n_; ++row ) {
					for( std::size_t column( 0 ); column < n_; ++column ) {
						isomorphism_address += connectivity( row, column ) * value_matrix( permuter[row], permuter[column] );
					}
				}
				map[ isomorphism_address ] = element_address;
			}
			++element_address;
		}
		
		// move to the next subgraph
		std::size_t row( 0 );
		std::size_t column( 1 );
		std::size_t * cell( &connectivity( row, column ) );
		++(*cell);
		while( *cell == r_pset ) {
			*cell = 0;
			if( d_ ) {
				std::swap( row, column );
				if( row < column ) {
					++column;
					if( column >= n_ ) {
						++row;
						column = row + 1;
						if( row == n_ - 1 ) {
							row = 0;
							column = 0;
						}
					}
				}
			} else {
				++column;
				if( column >= n_ ) {
					++row;
					column = row + 1;
					if( row == n_ - 1 ) {
						row = 0;
						column = 0;
					}
				}
			}
			cell = &connectivity( row, column );
			++(*cell);
		}
	}
}

template <std::size_t nn>
std::size_t vcp_static_mapper::subgraph_address( square_matrix<std::size_t,nn> const & connectivity ) const {
	std::size_t address( 0 );
	for( std::size_t row( 0 ); row < nn; ++row ) {
		for( std::size_t column( d_ ? 0 : row + 1 ); column < nn; ++column ) {
			address += connectivity( row, column ) * value_matrix( row, column );
		}
	}
	return address;
}

template <std::size_t nn>
std::size_t vcp_static_mapper::element_address( square_matrix<std::size_t,nn> const & connectivity ) const {
	return map[ subgraph_address( connectivity ) ];
}

std::size_t vcp_static_mapper::subgraph_address( square_matrix<std::size_t> const & connectivity ) const {
	std::size_t address( 0 );
	for( std::size_t row( 0 ); row < n_; ++row ) {
		for( std::size_t column( d_ ? 0 : row + 1 ); column < n_; ++column ) {
			address += connectivity( row, column ) * value_matrix( row, column );
		}
	}
	return address;
}

std::size_t vcp_static_mapper::element_address( square_matrix<std::size_t> const & connectivity ) const {
	return map[ subgraph_address( connectivity ) ];
}

square_matrix<std::size_t> vcp_static_mapper::element_structure( std::size_t address ) const {
	std::size_t r_pset( std::pow( 2, r_ ) );
	square_matrix<std::size_t> matrix( n_ );
	std::size_t i = 0;
	std::size_t j = 1;
	while( address > 0 ) {
		matrix( i, j ) = address % r_pset;
		address /= r_pset;
		if( d_ ) {
			std::swap( i, j );
			if( i < j ) {
				++j;
				if( j >= n_ ) {
					++i;
					j = i + 1;
				}
			}
		} else {
			++j;
			if( j >= n_ ) {
				++i;
				j = i + 1;
			}
		}
	}
	return matrix;
}

std::ostream & operator<<( std::ostream & os, vcp_static_mapper const & mapper ) {
	std::copy( mapper.map.begin(), mapper.map.end(), std::ostream_iterator<std::size_t>( os, "\n" ) );
	return os;
}

}

#endif

