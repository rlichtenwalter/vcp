/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_MULTIRELATIONAL_DIRECTED_GRAPH
#define VCP_MULTIRELATIONAL_DIRECTED_GRAPH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <vcp/multirelational_graph.hpp>

namespace vcp {

typedef std::size_t vertex_id_t;
typedef std::size_t edge_id_t;
typedef void * const * const_vertex_iterator;
typedef void * const * const_edge_iterator;

template <std::size_t r>
class multirelational_directed_graph {
	public:
		typedef typename multirelational_graph<r>::connectivity_address_type connectivity_address_type;
		multirelational_directed_graph();
		multirelational_directed_graph( multirelational_directed_graph const & );
		~multirelational_directed_graph();
		multirelational_directed_graph & operator=( multirelational_directed_graph const & );
		std::size_t vertex_count() const;
		std::size_t out_edge_count() const;
		std::size_t in_edge_count() const;
		std::size_t relation_count() const;
		const_vertex_iterator vertices_begin() const;
		const_vertex_iterator vertices_end() const;
		const_edge_iterator out_edges_begin() const;
		const_edge_iterator out_edges_end() const;
		const_edge_iterator in_edges_begin() const;
		const_edge_iterator in_edges_end() const;
		const_edge_iterator out_neighbors_begin( const_vertex_iterator ) const;
		const_edge_iterator out_neighbors_end( const_vertex_iterator ) const;
		const_edge_iterator in_neighbors_begin( const_vertex_iterator ) const;
		const_edge_iterator in_neighbors_end( const_vertex_iterator ) const;
		vertex_id_t vertex_id( const_vertex_iterator ) const;
		const_vertex_iterator target_of( const_edge_iterator ) const;
		edge_id_t edge_id( const_edge_iterator ) const;
		bool edge_exists( const_edge_iterator ) const;
		connectivity_address_type edge_value( const_edge_iterator ) const;
		const_edge_iterator out_edge( const_vertex_iterator, const_vertex_iterator ) const;
		const_edge_iterator in_edge( const_vertex_iterator, const_vertex_iterator ) const;
		bool out_edge_exists( const_vertex_iterator, const_vertex_iterator ) const;
		bool in_edge_exists( const_vertex_iterator, const_vertex_iterator ) const;
		template <std::size_t r_> friend std::ostream & operator<<( std::ostream &, multirelational_directed_graph<r_> const & );
		template <std::size_t r_> friend std::istream & operator>>( std::istream &, multirelational_directed_graph<r_> & );
	private:
		std::size_t num_vertices;
		std::size_t num_out_edges;
		std::unique_ptr<void*[]> vertices;
		std::unique_ptr<void*[]> edges;
		std::unique_ptr<connectivity_address_type[]> edge_values;
};

template <std::size_t r>
multirelational_directed_graph<r>::multirelational_directed_graph() : num_vertices(0), num_out_edges(0), vertices(std::unique_ptr<void*[]>(new void*[1])), edges(std::unique_ptr<void*[]>(new void*[1])), edge_values( std::unique_ptr<typename multirelational_directed_graph<r>::connectivity_address_type[]>(new typename multirelational_directed_graph<r>::connectivity_address_type[0])) {
	vertices[0] = &edges[0];
	edges[0] = NULL;
}

template <std::size_t r>
multirelational_directed_graph<r>::multirelational_directed_graph( multirelational_directed_graph const & g ) : num_vertices(g.num_vertices), num_out_edges(g.num_out_edges), vertices(std::unique_ptr<void*[]>(new void*[2*g.vertex_count()+1])), edges(std::unique_ptr<void*[]>(new void*[g.out_edge_count()+g.in_edge_count()+1])), edge_values( std::unique_ptr<typename multirelational_directed_graph<r>::connectivity_address_type[]>(new typename multirelational_directed_graph<r>::connectivity_address_type[2*g.num_out_edges])) {
	for( const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it ) {
		vertices[ g.vertex_id( it ) ] = &edges[ g.edge_id( g.out_neighbors_begin( it ) ) ];
		vertices[ vertex_count() + g.vertex_id( it ) ] = &edges[ g.edge_id( g.in_neighbors_begin( it ) ) ];
	}
	vertices[ 2*vertex_count() ]  = &edges[ out_edge_count() + in_edge_count() ];
	for( const_edge_iterator it = g.out_edges_begin(); it != g.out_edges_end(); ++it ) {
		edges[ g.edge_id( it ) ] = &vertices[ g.vertex_id( g.target_of( it ) ) ];
		edge_values[ g.edge_id( it ) ] = g.edge_value( it );
	}
	for( const_edge_iterator it = g.in_edges_begin(); it != g.in_edges_end(); ++it ) {
		edges[ g.edge_id( it ) ] = &vertices[ g.vertex_id( g.target_of( it ) ) ];
		edge_values[ g.edge_id( it ) ] = g.edge_value( it );
	}
	edges[ out_edge_count() + in_edge_count() ] = NULL;
}

template <std::size_t r>
multirelational_directed_graph<r>::~multirelational_directed_graph() {
}

template <std::size_t r>
multirelational_directed_graph<r> & multirelational_directed_graph<r>::operator=( multirelational_directed_graph const & g ) {
	if( this != &g ) {
		num_vertices = g.num_vertices;
		num_out_edges = g.num_out_edges;
		vertices = std::unique_ptr<void*[]>(new void*[ 2*g.vertex_count()+1 ]);
		edges = std::unique_ptr<void*[]>(new void*[ g.out_edge_count()+g.in_edge_count()+1 ]);
		edge_values = std::unique_ptr<typename multirelational_directed_graph<r>::connectivity_address_type[]>(new typename multirelational_directed_graph<r>::connectivity_address_type[2*g.num_out_edges]);
		for( const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it ) {
			vertices[ g.vertex_id( it ) ] = &edges[ g.edge_id( g.out_neighbors_begin( it ) ) ];
			vertices[ vertex_count() + g.vertex_id( it ) ] = &edges[ g.edge_id( g.in_neighbors_begin( it ) ) ];
		}
		vertices[ 2*vertex_count() ]  = &edges[ out_edge_count() + in_edge_count() ];
		for( const_edge_iterator it = g.out_edges_begin(); it != g.out_edges_end(); ++it ) {
			edges[ g.edge_id( it ) ] = &vertices[ g.vertex_id( g.target_of( it ) ) ];
			edge_values[ g.edge_id( it ) ] = g.edge_value( it );

		}
		for( const_edge_iterator it = g.in_edges_begin(); it != g.in_edges_end(); ++it ) {
			edges[ g.edge_id( it ) ] = &vertices[ g.vertex_id( g.target_of( it ) ) ];
			edge_values[ g.edge_id( it ) ] = g.edge_value( it );
		}
		edges[ out_edge_count() + in_edge_count() ] = NULL;
	}
	return *this;
}

template <std::size_t r>
std::size_t multirelational_directed_graph<r>::vertex_count() const {
	return num_vertices;
}

template <std::size_t r>
std::size_t multirelational_directed_graph<r>::out_edge_count() const {
	return num_out_edges;
}

template <std::size_t r>
std::size_t multirelational_directed_graph<r>::in_edge_count() const {
	return num_out_edges;
}

template <std::size_t r>
std::size_t multirelational_directed_graph<r>::relation_count() const {
	static std::size_t count( std::ceil( std::log2( *std::max_element( &edge_values[0], &edge_values[2*num_out_edges] ) + 1 ) ) );
	return count;
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::vertices_begin() const {
	return &vertices[ 0 ];
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::vertices_end() const {
	return &vertices[ vertex_count() ];
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::out_edges_begin() const {
	return &edges[ 0 ];
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::out_edges_end() const {
	return &edges[ out_edge_count() ];
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::in_edges_begin() const {
	return &edges[ out_edge_count() ];
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::in_edges_end() const {
	return &edges[ out_edge_count() + in_edge_count() ];
}

template <std::size_t r>
const_edge_iterator multirelational_directed_graph<r>::out_neighbors_begin( const_vertex_iterator it ) const {
	return static_cast<const_edge_iterator>( *it );
}

template <std::size_t r>
const_edge_iterator multirelational_directed_graph<r>::out_neighbors_end( const_vertex_iterator it ) const {
	return static_cast<const_edge_iterator>( *(it+1) );
}

template <std::size_t r>
const_edge_iterator multirelational_directed_graph<r>::in_neighbors_begin( const_vertex_iterator it ) const {
	return static_cast<const_edge_iterator>( *(vertex_count()+it) );
}

template <std::size_t r>
const_edge_iterator multirelational_directed_graph<r>::in_neighbors_end( const_vertex_iterator it ) const {
	return static_cast<const_edge_iterator>( *(vertex_count()+it+1) );
}

template <std::size_t r>
vertex_id_t multirelational_directed_graph<r>::vertex_id( const_vertex_iterator it ) const {
	return it - static_cast<const_vertex_iterator>( vertices.get() );
}

template <std::size_t r>
const_vertex_iterator multirelational_directed_graph<r>::target_of( const_edge_iterator it ) const {
	return static_cast<const_vertex_iterator>( *it );
}

template <std::size_t r>
edge_id_t multirelational_directed_graph<r>::edge_id( const_edge_iterator it ) const {
	return it - static_cast<const_edge_iterator>( edges.get() );
}

template <std::size_t r>
bool multirelational_directed_graph<r>::edge_exists( const_edge_iterator it ) const {
	return it != in_edges_end();
}

template <std::size_t r>
typename multirelational_directed_graph<r>::connectivity_address_type multirelational_directed_graph<r>::edge_value( const_edge_iterator it ) const {
	return edge_values[ edge_id( it ) ];
}

template <std::size_t r>
const_edge_iterator multirelational_directed_graph<r>::out_edge( const_vertex_iterator source, const_vertex_iterator target ) const {
	const_edge_iterator it( std::find( out_neighbors_begin( source ), out_neighbors_end( source ), target ) );
	return it == out_neighbors_end( source ) ? in_edges_end() : it;
}

template <std::size_t r>
const_edge_iterator multirelational_directed_graph<r>::in_edge( const_vertex_iterator source, const_vertex_iterator target ) const {
	const_edge_iterator it( std::find( in_neighbors_begin( source ), in_neighbors_end( source ), target ) );
	return it == in_neighbors_end( source ) ? in_edges_end() : it;
}

template <std::size_t r>
bool multirelational_directed_graph<r>::out_edge_exists( const_vertex_iterator source, const_vertex_iterator target ) const {
	return out_neighbors_end( source ) != std::find( out_neighbors_begin( source ), out_neighbors_end( source ), target );
}

template <std::size_t r>
bool multirelational_directed_graph<r>::in_edge_exists( const_vertex_iterator source, const_vertex_iterator target ) const {
	return in_neighbors_end( source ) != std::find( in_neighbors_begin( source ), in_neighbors_end( source ), target );
}

template <std::size_t r>
std::ostream & operator<<( std::ostream & os, multirelational_directed_graph<r> const & g ) {
	for( const_vertex_iterator vIt( g.vertices_begin() ); vIt < g.vertices_end(); ++vIt ) {
		const_edge_iterator nIt( g.out_neighbors_begin( vIt ) );
		while( nIt < g.out_neighbors_end( vIt ) - 1 ) {
			os << g.vertex_id( g.target_of( nIt ) ) << ',' << g.edge_value( nIt ) << ' ';
			++nIt;
		}
		if( nIt < g.out_neighbors_end( vIt ) ) {
			os << g.vertex_id( g.target_of( nIt ) ) << ',' << g.edge_value( nIt );
			++nIt;
		}
		if( vIt < g.vertices_end() ) {
			os << '\n';
		}
	}
	return os;
}

template <std::size_t r>
std::istream & operator>>( std::istream & is, multirelational_directed_graph<r> & g ) {
	std::vector<edge_id_t> out_v_temp;
	std::vector<std::pair<vertex_id_t,typename multirelational_directed_graph<r>::connectivity_address_type> > out_e_temp;
	std::string s1;
	std::size_t edge_id( 0 );
	
	while( std::getline( is, s1, '\n' ) ) {
		out_v_temp.push_back( edge_id );
		std::istringstream iss( s1 );
		std::string s2;
		while( std::getline( iss, s2, ' ' ) ) {
			std::istringstream iss2( s2 );
			std::string neighbor_str;
			std::string value_str;
			getline( iss2, neighbor_str, ',' );
			getline( iss2, value_str, ',' );
			std::istringstream iss3( neighbor_str );
			std::istringstream iss4( value_str );
			vertex_id_t neighbor;
			typename multirelational_directed_graph<r>::connectivity_address_type value;
			iss3 >> neighbor;
			iss4 >> value;
			out_e_temp.push_back( std::make_pair( neighbor, value ) );
			++edge_id;
		}
	}

	g.num_vertices = out_v_temp.size();
	g.num_out_edges = out_e_temp.size();
	
	g.vertices = std::unique_ptr<void*[]>(new void*[ 2 * g.vertex_count() + 1 ]);
	g.edges = std::unique_ptr<void*[]>(new void*[ g.out_edge_count() + g.in_edge_count() + 1 ]);
	g.edge_values = std::unique_ptr<typename multirelational_directed_graph<r>::connectivity_address_type[]>(new typename multirelational_directed_graph<r>::connectivity_address_type[g.out_edge_count()+g.in_edge_count()]);
	
	for( size_t i = 0; i < g.vertex_count(); ++i ) {
		g.vertices[ i ] = &g.edges[ out_v_temp[ i ] ];
	}
	g.vertices[ g.vertex_count() ] = &g.edges[ g.out_edge_count() ];
	for( size_t i( 0 ); i < g.out_edge_count(); ++i ) {
		g.edges[ i ] = &g.vertices[ out_e_temp[ i ].first ];
		g.edge_values[ i ] = out_e_temp[ i ].second;
	}

	out_v_temp = std::vector<edge_id_t>();
	out_e_temp = std::vector<std::pair<vertex_id_t,typename multirelational_directed_graph<r>::connectivity_address_type> >();

	std::vector<std::vector<std::pair<vertex_id_t,typename multirelational_directed_graph<r>::connectivity_address_type> > > in_temp( g.vertex_count() );
	for( const_vertex_iterator vIt( g.vertices_begin() ); vIt != g.vertices_end(); ++vIt ) {
		for( const_edge_iterator eIt( g.out_neighbors_begin( vIt ) ); eIt != g.out_neighbors_end( vIt ); ++eIt ) {
			in_temp[ g.vertex_id( g.target_of( eIt ) ) ].push_back( std::make_pair( g.vertex_id( vIt ), g.edge_value( eIt ) ) );
		}
	}

	std::size_t in_count( 0 );
	for( size_t i( 0 ); i < g.vertex_count(); ++i ) {
		g.vertices[ g.vertex_count()+i ] = &g.edges[ g.out_edge_count() + in_count ];
		std::vector<std::pair<vertex_id_t,typename multirelational_directed_graph<r>::connectivity_address_type> > const & in_neighbors = in_temp[ i ];
		for( typename std::vector<std::pair<vertex_id_t,typename multirelational_directed_graph<r>::connectivity_address_type> >::const_iterator it( in_neighbors.begin() ); it != in_neighbors.end(); ++it ) {
			g.edges[ g.out_edge_count() + in_count ] = &g.vertices[ it->first ];
			g.edge_values[ g.out_edge_count() + in_count ] = it->second;
			++in_count;
		}
	}

	g.vertices[ 2 * g.vertex_count() ] = &g.edges[ g.out_edge_count() + g.in_edge_count() ];
	g.edges[ g.out_edge_count() + g.in_edge_count() ] = NULL;
	
	return is;
}

}

#endif
