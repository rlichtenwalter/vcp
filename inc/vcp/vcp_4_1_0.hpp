/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_VCP_4_1_0
#define VCP_VCP_4_1_0

#include <array>
#include <cassert>
#include <cstddef>
#include <memory>
#include <utility>
#include <vcp/graph.hpp>

namespace vcp {

template <std::size_t n,std::size_t r,bool d> class vcp;

template <>
class vcp<4,1,0> {
	private:
		constexpr static const std::size_t num_elements = 40;
	public:
		vcp( graph const & g );
		constexpr static std::size_t element_count();
		std::array<unsigned long,num_elements> const generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 );
	private:
		enum connectivity_value {
			V1V2 = 1,
			V1V3 = 2,
			V1V4 = 4,
			V2V3 = 8,
			V2V4 = 16,
			V3V4 = 32
		};
		graph const & g;
		constexpr static const std::size_t num_structures = 64;
		static std::size_t element_address( std::size_t subgraph_address );
		unsigned long unconnected_pairs;
		std::unique_ptr<std::pair<const_vertex_iterator,unsigned char>[]> v3Vertices;
};

std::size_t vcp<4,1,0>::element_address( std::size_t subgraph_address ) {
	constexpr static const std::array<std::size_t,num_structures> map = {{0,1,2,3,2,3,4,5,6,7,8,9,10,11,12,13,6,7,10,11,8,9,12,13,14,15,16,17,16,17,18,19,20,21,22,23,22,23,24,25,26,27,28,29,30,31,32,33,26,27,30,31,28,29,32,33,34,35,36,37,36,37,38,39}};
	return map[ subgraph_address ];
}

constexpr std::size_t vcp<4,1,0>::element_count() {
	return num_elements;
}

vcp<4,1,0>::vcp( graph const & g ) :
		g( g ),
		unconnected_pairs( (g.vertex_count() * (g.vertex_count() - 1) / 2) - g.edge_count() ),
		v3Vertices( std::unique_ptr<std::pair<const_vertex_iterator,unsigned char>[]>(new std::pair<const_vertex_iterator,unsigned char>[ MAX_NEIGHBORS ] )) {
}


std::array<unsigned long,vcp<4,1,0>::element_count()> const vcp<4,1,0>::generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 ) {
	std::array<unsigned long,element_count()> counts = {{0}};
	
	std::size_t v1v2( V1V2 * g.edge_exists( v1, v2 ) );
	
	unsigned long connections( 0 );
	unsigned long gaps( 0 );

	const_edge_iterator v1_neighbors_it( g.neighbors_begin( v1 ) );
	const_edge_iterator v1_neighbors_end( g.neighbors_end( v1 ) );
	const_edge_iterator v2_neighbors_it( g.neighbors_begin( v2 ) );
	const_edge_iterator v2_neighbors_end( g.neighbors_end( v2 ) );
	assert( MAX_NEIGHBORS > (v1_neighbors_end - v1_neighbors_it) + (v2_neighbors_end - v2_neighbors_it ) ); // this should always be contiguous storage; we can only over-allocate by a factor of 2, which is of much lower cost than maintaining a doubly-linked list; there exists a strict upper bound on the final size of v3Vertices
	std::pair<const_vertex_iterator,unsigned char>* v3Vertices_begin( &v3Vertices[0] );
	std::pair<const_vertex_iterator,unsigned char>* v3Vertices_end( &v3Vertices[0] );
	while( v1_neighbors_it != v1_neighbors_end && v2_neighbors_it != v2_neighbors_end ) {
		if( g.target_of( v1_neighbors_it ) < g.target_of( v2_neighbors_it )  ) {
			if( g.target_of( v1_neighbors_it ) != v2 ) {
				++connections;
				++gaps;
				v3Vertices_end->first = g.target_of( v1_neighbors_it );
				v3Vertices_end->second = v1v2 + V1V3;
				++v3Vertices_end;
			}
			++v1_neighbors_it;
		} else if( g.target_of( v1_neighbors_it ) > g.target_of( v2_neighbors_it )  ) {
			if( g.target_of( v2_neighbors_it ) != v1 ) {
				++connections;
				++gaps;
				v3Vertices_end->first = g.target_of( v2_neighbors_it );
				v3Vertices_end->second = v1v2 + V2V3;
				++v3Vertices_end;
			}
			++v2_neighbors_it;
		} else { // the next neighbor is shared by both v1 and v2, so it cannot be either and we do not need to check to exclude it
			connections += 2;
			v3Vertices_end->first = g.target_of( v1_neighbors_it );
			v3Vertices_end->second = v1v2 + V1V3 + V2V3;
			++v3Vertices_end;
			++v1_neighbors_it;
			++v2_neighbors_it;
		}
	}
	while( v1_neighbors_it != v1_neighbors_end ) {
		if( g.target_of( v1_neighbors_it ) != v2 ) {
			++connections;
			++gaps;
			v3Vertices_end->first = g.target_of( v1_neighbors_it );
			v3Vertices_end->second = v1v2 + V1V3;
			++v3Vertices_end;
		}
		++v1_neighbors_it;
	}
	while( v2_neighbors_it != v2_neighbors_end ) {
		if( g.target_of( v2_neighbors_it ) != v1 ) {
			++connections;
			++gaps;
			v3Vertices_end->first = g.target_of( v2_neighbors_it );
			v3Vertices_end->second = v1v2 + V2V3;
			++v3Vertices_end;
		}
		++v2_neighbors_it;
	}
	
	unsigned long v3_count( v3Vertices_end - v3Vertices_begin );
	unsigned long v4_count( 0 );
	for( std::pair<const_vertex_iterator,unsigned char>* it1( v3Vertices_begin ); it1 != v3Vertices_end; ++it1 ) { // for each v3 vertex computed above
		const_edge_iterator v3_neighbors_it( g.neighbors_begin( it1->first ) );
		const_edge_iterator v3_neighbors_end( g.neighbors_end( it1->first ) );
		unsigned long v4_local_count( 0 ); // keep track of how many v4 vertices are only the result of the neighbors of this v3
		for( std::pair<const_vertex_iterator,unsigned char>* it2( v3Vertices_begin ); it2 != v3Vertices_end; ++it2 ) { // consider other v3 vertices as candidate v4 vertices
			while( v3_neighbors_it != v3_neighbors_end && g.target_of( v3_neighbors_it ) < it2->first ) { // the v3 neighbor is exclusively a v4 vertex
				if( g.target_of( v3_neighbors_it ) != v1 && g.target_of( v3_neighbors_it ) != v2 ) { // if this exclusively v4 vertex is not v1 or v2
					++v4_local_count;
				}
				++v3_neighbors_it;
			}
			if( v3_neighbors_it == v3_neighbors_end || g.target_of( v3_neighbors_it ) > it2->first ) { // there is no edge between the v3 vertex and the other v3 vertex serving as a v4 vertex
				if( it1->first < it2->first ) { // to be a candidate vertex, the other v3 vertex must be greater to avoid double counting
					++gaps;
					++counts[ element_address( it1->second + (it2->second == v1v2 + V1V3 ? V1V4 : (it2->second == v1v2 + V2V3 ? V2V4 : V1V4 + V2V4)) ) ];
				}
			} else { // there is an edge between the v3 vertex and the other v3 vertex serving as a v4 vertex
				if( it1->first < it2->first ) { // to be a candidate vertex, the other v3 vertex must be greater to avoid double counting
					++connections;
					++counts[ element_address( it1->second + (it2->second == v1v2 + V1V3 ? V1V4 : (it2->second == v1v2 + V2V3 ? V2V4 : V1V4 + V2V4)) + V3V4 ) ];
				}
				++v3_neighbors_it;
			}
		}
		while( v3_neighbors_it != v3_neighbors_end ) { // we have to be sure to go through the rest of the neighbors of v3, and all of these are exclusively v4
			if( g.target_of( v3_neighbors_it ) != v1 && g.target_of( v3_neighbors_it ) != v2 ) {
				++v4_local_count;
			}
			++v3_neighbors_it;
		}
		v4_count += v4_local_count;
		connections += v4_local_count;
		gaps += 2*v4_local_count;
		counts[ element_address( it1->second + V3V4 ) ] += v4_local_count;
		counts[ element_address( it1->second ) ] += g.vertex_count() - 2 - v3_count - v4_local_count;
	}

	// account for the least connected substructures
	counts[ element_address( v1v2+V3V4 ) ] = g.edge_count() - (connections + static_cast<bool>(v1v2));
	
	// the final computation is somewhat complicated
	// we start with the number of total gaps in the network, subtract the number of directly observed gaps, then subtract the number of gaps contributed by triangles that we do not directly observe that result from vertices that we do not directly observe
	// the number of unobserved triangles is recorded above as g.vertex_count() - 2 - v3_count - v4_local_count
	// outside the loop this is expressed as v3_count * (g.vertex_count() - 2 - v3_count) - v4_count
	// this is also the number of unobserved gaps that result from the v3 vertices in combination with vertices we do not directly observe
	// but we also have to account for the unobserved gaps between v1 and v2 and these unobserved vertices
	// this is pretty clearly 2*(g.vertex_count() - 2 - v3_count - v4_count), a gap for each pairing of unobserved vertex with v1 or v2
	// the computation below is thus equivalent to:
	// unconnected_pairs - (gaps + !static_cast<bool>(v1v2)) - 2*(g.vertex_count() - 2 - v3_count - v4_count) - v3_count*(g.vertex_count() - 2 - v3_count) + v4_count
	// we can simplify this expression as below
	counts[ element_address(v1v2) ] = unconnected_pairs - (gaps + !static_cast<bool>(v1v2)) - (2 + v3_count) * (g.vertex_count() - 2 - v3_count) + 3 * v4_count;
	
	return counts;
}

}

#endif
