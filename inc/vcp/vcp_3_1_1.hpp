/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_VCP_3_1_1
#define VCP_VCP_3_1_1

#include <array>
#include <cstddef>
#include <vcp/directed_graph.hpp>

namespace vcp {

template <std::size_t n,std::size_t r,bool d> class vcp;

template <>
class vcp<3,1,1> {
	private:
		const static unsigned long num_elements = 64;
	public:
		vcp( directed_graph const & g );
		constexpr static std::size_t element_count();
		std::array<unsigned long,num_elements> const generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 );
	private:
		enum directedness_value {
			OUT = 1,
			IN = 2,
			BOTH = 3
		};
		enum connectivity_value {
			V1V2 = 1,
			V1V3 = 4,
			V2V3 = 16
		};
		directed_graph const & g;
		std::pair<const_edge_iterator,directedness_value> next_union_element( const_edge_iterator & it1, const_edge_iterator end1, const_edge_iterator & it2, const_edge_iterator end2 ) const;
};

constexpr std::size_t vcp<3,1,1>::element_count() {
	return num_elements;
}

vcp<3,1,1>::vcp( directed_graph const & g ) : g(g) {
}

std::pair<const_edge_iterator,vcp<3,1,1>::directedness_value> vcp<3,1,1>::next_union_element( const_edge_iterator & it1, const_edge_iterator end1, const_edge_iterator & it2, const_edge_iterator end2 ) const { // out-neighbor iterators should always come first
	if( it1 != end1 && it2 != end2 ) {
		if( g.target_of( it1 ) < g.target_of( it2 ) ) {
			return std::make_pair( it1++, OUT );
		} else if( g.target_of( it1 ) > g.target_of( it2 ) ) {
			return std::make_pair( it2++, IN );
		} else {
			++it2;
			return std::make_pair( it1++, BOTH );
		}
	} else if( it1 != end1 ) {
			return std::make_pair( it1++, OUT );
	} else if( it2 != end2 ) {
			return std::make_pair( it2++, IN );
	} else {
		return std::make_pair( end2, OUT ); // the second element of this pair should never be used
	}
}

std::array<unsigned long,vcp<3,1,1>::element_count()> const vcp<3,1,1>::generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 ) {
	std::array<unsigned long,element_count()> counts = {{0}};
	
	std::size_t v1v2( g.out_edge_exists( v1, v2 ) * OUT + g.in_edge_exists( v1, v2 ) * IN );
	
	const_edge_iterator v1_out_neighbors_it( g.out_neighbors_begin( v1 ) );
	const_edge_iterator v1_out_neighbors_end( g.out_neighbors_end( v1 ) );
	const_edge_iterator v1_in_neighbors_it( g.in_neighbors_begin( v1 ) );
	const_edge_iterator v1_in_neighbors_end( g.in_neighbors_end( v1 ) );
	const_edge_iterator v2_out_neighbors_it( g.out_neighbors_begin( v2 ) );
	const_edge_iterator v2_out_neighbors_end( g.out_neighbors_end( v2 ) );
	const_edge_iterator v2_in_neighbors_it( g.in_neighbors_begin( v2 ) );
	const_edge_iterator v2_in_neighbors_end( g.in_neighbors_end( v2 ) );
	std::pair<const_edge_iterator,directedness_value> min1( next_union_element( v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it, v1_in_neighbors_end ) );
	std::pair<const_edge_iterator,directedness_value> min2( next_union_element( v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it, v2_in_neighbors_end ) );
	unsigned long union_cardinality( 0 );
	while( min1.first != v1_in_neighbors_end && min2.first != v2_in_neighbors_end ) {
		if( g.target_of( min1.first ) < g.target_of( min2.first ) ) {
			if( g.target_of( min1.first ) != v2 ) {
				++union_cardinality;
				++counts[ v1v2 + min1.second * V1V3 ];
			}
			min1 = next_union_element( v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it, v1_in_neighbors_end );
		} else if( g.target_of( min1.first ) > g.target_of( min2.first ) ) {
			if( g.target_of( min2.first ) != v1 ) {
				++union_cardinality;
				++counts[ v1v2 + min2.second * V2V3 ];
			}
			min2 = next_union_element( v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it, v2_in_neighbors_end );
		} else { // the next neighbor is shared by both v1 and v2, so it cannot be either and we do not need to check to exclude it
			++union_cardinality;
			++counts[ v1v2 + min1.second * V1V3 + min2.second * V2V3 ];
			min1 = next_union_element( v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it, v1_in_neighbors_end );
			min2 = next_union_element( v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it, v2_in_neighbors_end );
		}
	}
	while( min1.first != v1_in_neighbors_end ) {
		if( g.target_of( min1.first ) != v2 ) {
			++union_cardinality;
			++counts[ v1v2 + min1.second * V1V3 ];
		}
		min1 = next_union_element( v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it, v1_in_neighbors_end );
	}
	while( min2.first != v2_in_neighbors_end ) {
		if( g.target_of( min2.first ) != v1 ) {
			++union_cardinality;
			++counts[ v1v2 + min2.second * V2V3 ];
		}
		min2 = next_union_element( v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it, v2_in_neighbors_end );
	}
	
	counts[ v1v2 ] = g.vertex_count() - 2 - union_cardinality;

	return counts;
}

}

#endif
