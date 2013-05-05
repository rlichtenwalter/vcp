/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_VCP
#define VCP_VCP

#include <cstddef>
#include <iostream>
#include <map>
#include <type_traits>
#include <utility>
#include <vcp/directed_graph.hpp>
#include <vcp/graph.hpp>
#include <vcp/multirelational_directed_graph.hpp>
#include <vcp/multirelational_graph.hpp>
#include <vcp/square_matrix.hpp>
#include <vcp/vcp_3_1_0.hpp>
#include <vcp/vcp_3_r_0.hpp>
#include <vcp/vcp_3_1_1.hpp>
#include <vcp/vcp_3_r_1.hpp>
#include <vcp/vcp_4_1_0.hpp>
#include <vcp/vcp_4_1_1.hpp>
#include <vcp/vcp_4_r_0.hpp>
#include <vcp/vcp_4_r_1.hpp>
#include <vcp/vcp_dynamic_mapper.hpp>

namespace vcp {
	
template <std::size_t n,std::size_t r,bool d>
class vcp {
	public:
		typedef typename std::conditional< d,
				typename std::conditional< (r>1), multirelational_directed_graph<r>, directed_graph >::type,
				typename std::conditional< (r>1), multirelational_graph<r>, graph >::type
				>::type graph_type;
		typedef typename multirelational_graph<r>::connectivity_address_type connectivity_address_type;
		typedef typename vcp_dynamic_mapper<n,r,d>::subgraph_address_type subgraph_address_type;
		vcp( graph_type const & g );
		std::map<subgraph_address_type,unsigned long> const generate_vector( const_vertex_iterator, const_vertex_iterator );
	private:
		graph_type const & g;
		vcp_dynamic_mapper<n,r,d> mapper;
		void helper( std::array<const_vertex_iterator,n> & vertices, std::size_t current_index, square_matrix<connectivity_address_type,n> & connectivity, std::map<subgraph_address_type,unsigned long> & counts );
};

const_edge_iterator edge( graph const & g, const_vertex_iterator v1, const_vertex_iterator v2 );
const_edge_iterator edge( directed_graph const & g, const_vertex_iterator v1, const_vertex_iterator v2 );
template <std::size_t r> const_edge_iterator edge( multirelational_graph<r> const & g, const_vertex_iterator v1, const_vertex_iterator v2 );
template <std::size_t r> const_edge_iterator edge( multirelational_directed_graph<r> const & g, const_vertex_iterator v1, const_vertex_iterator v2 );
bool edge_value( graph const & g, const_edge_iterator edge );
bool edge_value( directed_graph const & g, const_edge_iterator edge );
template <std::size_t r> typename multirelational_graph<r>::connectivity_address_type edge_value( multirelational_graph<r> const & g, const_edge_iterator edge );
template <std::size_t r> typename multirelational_graph<r>::connectivity_address_type edge_value( multirelational_directed_graph<r> const & g, const_edge_iterator edge );

const_edge_iterator edge( graph const & g, const_vertex_iterator v1, const_vertex_iterator v2 ) {
	return g.edge( v1, v2 );
}

const_edge_iterator edge( directed_graph const & g, const_vertex_iterator v1, const_vertex_iterator v2 ) {
	return g.out_edge( v1, v2 );
}

template <std::size_t r>
const_edge_iterator edge( multirelational_graph<r> const & g, const_vertex_iterator v1, const_vertex_iterator v2 ) {
	return g.edge( v1, v2 );
}
template <std::size_t r>
const_edge_iterator edge( multirelational_directed_graph<r> const & g, const_vertex_iterator v1, const_vertex_iterator v2 ) {
	return g.out_edge( v1, v2 );
}

bool edge_value( graph const & g, const_edge_iterator edge ) {
	return g.edge_exists( edge );
}

bool edge_value( directed_graph const & g, const_edge_iterator edge ) {
	return g.edge_exists( edge );
}

template <typename std::size_t r>
typename multirelational_graph<r>::connectivity_address_type edge_value( multirelational_graph<r> const & g, const_edge_iterator edge ) {
	return g.edge_value( edge );
}

template <typename std::size_t r>
typename multirelational_graph<r>::connectivity_address_type edge_value( multirelational_directed_graph<r> const & g, const_edge_iterator edge ) {
	return g.edge_value( edge );
}

template <std::size_t n,std::size_t r,bool d>
vcp<n,r,d>::vcp( graph_type const & g ) : g( g ) {
}

template <std::size_t n,std::size_t r,bool d>
void vcp<n,r,d>::helper( std::array<const_vertex_iterator,n> & vertices, std::size_t current_index, square_matrix<connectivity_address_type,n> & connectivity, std::map<subgraph_address_type,unsigned long> & counts ) {
	if( vertices[ current_index ] == vertices[ 0 ] || vertices[ current_index ] == vertices[ 1 ] ) {
		return;
	}
	for( std::size_t index( 0 ); index < current_index; ++index ) {
		connectivity( index, current_index ) = edge_value( g, edge( g, vertices[ index ], vertices[ current_index ] ) );
		if( d ) {
			connectivity( current_index, index ) = edge_value( g, edge( g, vertices[ current_index ], vertices[ index ] ) );
		}
	}
	if( n == current_index + 1 ) {
		++counts.insert( std::make_pair( mapper.canonical_subgraph_address( connectivity ), 0 ) ).first->second;
	} else {
		for( const_vertex_iterator v( vertices[ current_index ] + 1 ); v != g.vertices_end(); ++v ) {
			vertices[ current_index + 1 ] = v;
			helper( vertices, current_index + 1, connectivity, counts );
		}
	}
	return;
}

template <std::size_t n,std::size_t r,bool d>
std::map<typename vcp<n,r,d>::subgraph_address_type,unsigned long> const vcp<n,r,d>::generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 ) {
	std::map<subgraph_address_type,unsigned long> counts;
	std::array<const_vertex_iterator,n> vertices;
	square_matrix<connectivity_address_type,n> connectivity;
	connectivity( 0, 1 ) = edge_value( g, edge( g, v1, v2 ) );
	if( d ) {
		connectivity( 1, 0 ) = edge_value( g, edge( g, v2, v1 ) );
	}
	vertices[ 0 ] = v1;
	vertices[ 1 ] = v2;
	for( const_vertex_iterator v3( g.vertices_begin() ); v3 != g.vertices_end(); ++v3 ) {
		vertices[ 2 ] = v3;
		helper( vertices, 2, connectivity, counts );
	}
	return counts;
}

}

#endif
