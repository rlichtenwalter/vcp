/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_VCP_4_R_0
#define VCP_VCP_4_R_0

#include <cstddef>
#include <cassert>
#include <map>
#include <utility>
#include <vcp/multirelational_graph.hpp>
#include <vcp/square_matrix.hpp>
#include <vcp/vcp_dynamic_mapper.hpp>

namespace vcp {

template <std::size_t n,std::size_t r,bool d> class vcp;

template <std::size_t r>
class vcp<4,r,0> {
	public:
		typedef typename multirelational_graph<r>::connectivity_address_type connectivity_address_type;
		typedef typename vcp_dynamic_mapper<4,r,0>::subgraph_address_type subgraph_address_type;
		vcp( multirelational_graph<r> const & g );
		std::map<subgraph_address_type,unsigned long> const generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 );
	private:
		typedef square_matrix<connectivity_address_type,4> connectivity_matrix;
		multirelational_graph<r> const & g;
		vcp_dynamic_mapper<4,r,0> mapper;
		std::map<connectivity_address_type,unsigned long> edge_types;
		std::unique_ptr<std::pair<const_vertex_iterator,connectivity_matrix>[]> v3Vertices;
};

template <std::size_t r>
vcp<4,r,0>::vcp( multirelational_graph<r> const & g ) : g( g ), mapper(), v3Vertices( std::unique_ptr<std::pair<const_vertex_iterator,connectivity_matrix>[]>(new std::pair<const_vertex_iterator,connectivity_matrix>[ MAX_NEIGHBORS ] )) {
	unsigned long & gaps( edge_types.insert( std::make_pair( 0, g.vertex_count() * (g.vertex_count() - 1) / 2 ) ).first->second );
	for( const_vertex_iterator it( g.vertices_begin() ); it != g.vertices_end(); ++it ) {
		for( const_edge_iterator eIt( g.neighbors_begin( it ) ); eIt != g.neighbors_end( it ); ++eIt ) {
			if( it < g.target_of( eIt ) ) {
				++edge_types.insert( std::make_pair( g.edge_value( eIt ), 0 ) ).first->second;
				--gaps;
			}
		}
	}
}

template <std::size_t r>
std::map<typename vcp<4,r,0>::subgraph_address_type,unsigned long> const vcp<4,r,0>::generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 ) {
	std::map<subgraph_address_type,unsigned long> counts;
	std::map<connectivity_address_type,unsigned long> temp_edge_types;

	connectivity_matrix connectivity;
	connectivity( 0, 1 ) = g.edge_value( g.edge( v1, v2 ) );

	unsigned long & gaps( temp_edge_types.insert( std::make_pair( 0, 0 ) ).first->second );

	const_edge_iterator v1_neighbors_it( g.neighbors_begin( v1 ) );
	const_edge_iterator v1_neighbors_end( g.neighbors_end( v1 ) );
	const_edge_iterator v2_neighbors_it( g.neighbors_begin( v2 ) );
	const_edge_iterator v2_neighbors_end( g.neighbors_end( v2 ) );
	assert( MAX_NEIGHBORS > (v1_neighbors_end - v1_neighbors_it) + (v2_neighbors_end - v2_neighbors_it ) ); // this should always be contiguous storage; we can only over-allocate by a factor of 2, which is of much lower cost than maintaining a doubly-linked list; there exists a strict upper bound on the final size of v3Vertices
	std::pair<const_vertex_iterator,connectivity_matrix>* v3Vertices_begin( &v3Vertices[0] );
	std::pair<const_vertex_iterator,connectivity_matrix>* v3Vertices_end( &v3Vertices[0] );
	while( v1_neighbors_it != v1_neighbors_end && v2_neighbors_it != v2_neighbors_end ) {
		if( g.target_of( v1_neighbors_it ) < g.target_of( v2_neighbors_it )  ) {
			if( g.target_of( v1_neighbors_it ) != v2 ) {
				++temp_edge_types.insert( std::make_pair( g.edge_value( v1_neighbors_it ), 0 ) ).first->second;
				++gaps;
				v3Vertices_end->first = g.target_of( v1_neighbors_it );
				v3Vertices_end->second = connectivity;
				v3Vertices_end->second( 0, 2 ) = g.edge_value( v1_neighbors_it );
				++v3Vertices_end;
			}
			++v1_neighbors_it;
		} else if( g.target_of( v1_neighbors_it ) > g.target_of( v2_neighbors_it )  ) {
			if( g.target_of( v2_neighbors_it ) != v1 ) {
				++temp_edge_types.insert( std::make_pair( g.edge_value( v2_neighbors_it ), 0 ) ).first->second;
				++gaps;
				v3Vertices_end->first = g.target_of( v2_neighbors_it );
				v3Vertices_end->second = connectivity;
				v3Vertices_end->second( 1, 2 ) = g.edge_value( v1_neighbors_it );
				++v3Vertices_end;
			}
			++v2_neighbors_it;
		} else { // the next neighbor is shared by both v1 and v2, so it cannot be either and we do not need to check to exclude it
			++temp_edge_types.insert( std::make_pair( g.edge_value( v1_neighbors_it ), 0 ) ).first->second;
			++temp_edge_types.insert( std::make_pair( g.edge_value( v2_neighbors_it ), 0 ) ).first->second;
			v3Vertices_end->first = g.target_of( v1_neighbors_it );
			v3Vertices_end->second = connectivity;
			v3Vertices_end->second( 0, 2 ) = g.edge_value( v1_neighbors_it );
			v3Vertices_end->second( 1, 2 ) = g.edge_value( v2_neighbors_it );
			++v3Vertices_end;
			++v1_neighbors_it;
			++v2_neighbors_it;
		}
	}
	while( v1_neighbors_it != v1_neighbors_end ) {
		if( g.target_of( v1_neighbors_it ) != v2 ) {
			++temp_edge_types.insert( std::make_pair( g.edge_value( v1_neighbors_it ), 0 ) ).first->second;
			++gaps;
			v3Vertices_end->first = g.target_of( v1_neighbors_it );
			v3Vertices_end->second = connectivity;
			v3Vertices_end->second( 0, 2 ) = g.edge_value( v1_neighbors_it );
			++v3Vertices_end;
		}
		++v1_neighbors_it;
	}
	while( v2_neighbors_it != v2_neighbors_end ) {
		if( g.target_of( v2_neighbors_it ) != v1 ) {
			++temp_edge_types.insert( std::make_pair( g.edge_value( v2_neighbors_it ), 0 ) ).first->second;
			++gaps;
			v3Vertices_end->first = g.target_of( v2_neighbors_it );
			v3Vertices_end->second = connectivity;
			v3Vertices_end->second( 1, 2 ) = g.edge_value( v1_neighbors_it );
			++v3Vertices_end;
		}
		++v2_neighbors_it;
	}
	
	std::size_t v3_count( v3Vertices_end-v3Vertices_begin );
	std::size_t v4_count( 0 );
	for( std::pair<const_vertex_iterator,connectivity_matrix>* it1( v3Vertices_begin ); it1 != v3Vertices_end; ++it1 ) { // for each v3 vertex computed above
		const_edge_iterator v3_neighbors_it( g.neighbors_begin( it1->first ) );
		const_edge_iterator v3_neighbors_end( g.neighbors_end( it1->first ) );
		unsigned long v4_local_count( 0 ); // keep track of how many v4 vertices are only the result of the neighbors of this v3
		for( std::pair<const_vertex_iterator,connectivity_matrix>* it2( v3Vertices_begin ); it2 != v3Vertices_end; ++it2 ) { // consider other v3 vertices as candidate v4 vertices
			while( v3_neighbors_it != v3_neighbors_end && g.target_of( v3_neighbors_it ) < it2->first ) { // the v3 neighbor is exclusively a v4 vertex
				if( g.target_of( v3_neighbors_it ) != v1 && g.target_of( v3_neighbors_it ) != v2 ) { // if this exclusively v4 vertex is not v1 or v2
					++temp_edge_types.insert( std::make_pair( g.edge_value( v3_neighbors_it ), 0 ) ).first->second;
					++v4_local_count;
					it1->second( 0, 3 ) = 0;
					it1->second( 1, 3 ) = 0;
					it1->second( 2, 3 ) = g.edge_value( v3_neighbors_it );
					++counts.insert( std::make_pair( mapper.canonical_subgraph_address( it1->second ), 0 ) ).first->second;
				}
				++v3_neighbors_it;
			}
			if( v3_neighbors_it == v3_neighbors_end || g.target_of( v3_neighbors_it ) > it2->first ) { // there is no edge between the v3 vertex and the other v3 vertex serving as a v4 vertex
				if( it1->first < it2->first ) { // to be a candidate vertex, the other v3 vertex must be greater to avoid double counting
					++gaps;
					it1->second( 0, 3 ) = it2->second( 0, 2 );
					it1->second( 1, 3 ) = it2->second( 1, 2 );
					it1->second( 2, 3 ) = 0;
					++counts.insert( std::make_pair( mapper.canonical_subgraph_address( it1->second ), 0 ) ).first->second;
				}
			} else { // there is an edge between the v3 vertex and the other v3 vertex serving as a v4 vertex
				if( it1->first < it2->first ) { // to be a candidate vertex, the other v3 vertex must be greater to avoid double counting
					++temp_edge_types.insert( std::make_pair( g.edge_value( v3_neighbors_it ), 0 ) ).first->second;
					it1->second( 0, 3 ) = it2->second( 0, 2 );
					it1->second( 1, 3 ) = it2->second( 1, 2 );
					it1->second( 2, 3 ) = g.edge_value( v3_neighbors_it );
					++counts.insert( std::make_pair( mapper.canonical_subgraph_address( it1->second ), 0 ) ).first->second;
				}
				++v3_neighbors_it;
			}
		}
		while( v3_neighbors_it != v3_neighbors_end ) { // we have to be sure to go through the rest of the neighbors of v3, and all of these are exclusively v4
			if( g.target_of( v3_neighbors_it ) != v1 && g.target_of( v3_neighbors_it ) != v2 ) {
				++temp_edge_types.insert( std::make_pair( g.edge_value( v3_neighbors_it ), 0 ) ).first->second;
				++v4_local_count;
				it1->second( 0, 3 ) = 0;
				it1->second( 1, 3 ) = 0;
				it1->second( 2, 3 ) = g.edge_value( v3_neighbors_it );
				++counts.insert( std::make_pair( mapper.canonical_subgraph_address( it1->second ), 0 ) ).first->second;
			}
			++v3_neighbors_it;
		}
		v4_count += v4_local_count;
		gaps += 2*v4_local_count;
		it1->second( 0, 3 ) = 0;
		it1->second( 1, 3 ) = 0;
		it1->second( 2, 3 ) = 0;
		counts.insert( std::make_pair( mapper.canonical_subgraph_address( it1->second ), 0 ) ).first->second += g.vertex_count() - 2 - v3_count - v4_local_count;
	}

	// account for the least connected substructures
	for( typename std::map<connectivity_address_type,unsigned long>::const_iterator it( edge_types.begin() ); it != edge_types.end(); ++it ) {
		connectivity( 2, 3 ) = it->first;
		unsigned long count = it->second;
		typename std::map<connectivity_address_type,unsigned long>::const_iterator temp_it( temp_edge_types.find( it->first ) );
		if( temp_it != temp_edge_types.end() ) {
			count -= temp_it->second;
			if( it->first == 0 ) {
				count -= !static_cast<bool>( connectivity( 0, 1 ) ) + (2 + v3_count) * (g.vertex_count() - 2 - v3_count) - 3 * v4_count;
			}
		}
		counts.insert( std::make_pair( mapper.canonical_subgraph_address( connectivity ), 0 ) ).first->second += count;
	}
	
	return counts;
}

}

#endif
