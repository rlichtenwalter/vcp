/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_VCP_3_1_0
#define VCP_VCP_3_1_0

#include <array>
#include <cstddef>
#include <vcp/graph.hpp>

namespace vcp {
	
template <std::size_t n,std::size_t r,bool d> class vcp;

template <>
class vcp<3,1,0> {
	private:
		constexpr static const std::size_t num_elements = 8;
	public:
		vcp( graph const & );
		constexpr static std::size_t element_count();
		std::array<unsigned long,num_elements> const generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 );
	private:
		enum connectivity_value {
			V1V2 = 1,
			V1V3 = 2,
			V2V3 = 4
		};
		graph const & g;
};

vcp<3,1,0>::vcp( graph const & g ) : g(g) {
}

constexpr std::size_t vcp<3,1,0>::element_count() {
	return num_elements;
}

std::array<unsigned long,vcp<3,1,0>::element_count()> const vcp<3,1,0>::generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 ) {
	std::array<unsigned long,element_count()> counts = {{0}};

	size_t v1v2( V1V2 * g.edge_exists( v1, v2 ) );

	counts[ v1v2 ] = g.vertex_count() - 2;
	const_edge_iterator v1_it( g.neighbors_begin( v1 ) );
	const_edge_iterator v1_end( g.neighbors_end( v1 ) );
	const_edge_iterator v2_it( g.neighbors_begin( v2 ) );
	const_edge_iterator v2_end( g.neighbors_end( v2 ) );
	
	while( v1_it != v1_end && v2_it != v2_end ) {
		if( g.target_of( v1_it ) == v2 ) {
			++v1_it;
		} else if( g.target_of( v2_it ) == v1 ) {
			++v2_it;
		} else {
			--counts[ v1v2 ];
			if( g.target_of( v1_it ) < g.target_of( v2_it ) ) {
				++counts[ v1v2 + V1V3 ];
				++v1_it;
			} else if( g.target_of( v1_it ) > g.target_of( v2_it ) ) {
				++counts[ v1v2 + V2V3 ];
				++v2_it;
			} else {
				++counts[ v1v2 + V1V3 + V2V3 ];
				++v1_it;
				++v2_it;
			}
		}
	}
	while( v1_it != v1_end ) {
		if( g.target_of( v1_it ) != v2 ) {
			--counts[ v1v2 ];
			++counts[ v1v2 + V1V3 ];
		}
		++v1_it;
	} while( v2_it != v2_end ) {
		if( g.target_of( v2_it ) != v1 ) {
			--counts[ v1v2 ];
			++counts[ v1v2 + V2V3 ];
		}
		++v2_it;
	}

	return counts;
}

}

#endif
