/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_MULTIRELATIONAL_GRAPH
#define VCP_MULTIRELATIONAL_GRAPH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
#include <vcp/graph.hpp>

namespace vcp {
	
typedef std::size_t vertex_id_t;
typedef std::size_t edge_id_t;
typedef void * const * const_vertex_iterator;
typedef void * const * const_edge_iterator;

template <std::size_t r>
class multirelational_graph {
	public:
		typedef typename std::conditional<r<=CHAR_BIT*sizeof(std::size_t),std::size_t,boost::multiprecision::number<boost::multiprecision::cpp_int_backend<r,r,boost::multiprecision::unsigned_magnitude,boost::multiprecision::unchecked,void> > >::type connectivity_address_type;
		multirelational_graph();
		multirelational_graph( multirelational_graph const & );
		~multirelational_graph();
		multirelational_graph & operator=( multirelational_graph const & );
		std::size_t vertex_count() const;
		std::size_t edge_count() const;
		std::size_t relation_count() const;
		const_vertex_iterator vertices_begin() const;
		const_vertex_iterator vertices_end() const;
		const_edge_iterator edges_begin() const;
		const_edge_iterator edges_end() const;
		const_edge_iterator neighbors_begin( const_vertex_iterator ) const;
		const_edge_iterator neighbors_end( const_vertex_iterator ) const;
		vertex_id_t vertex_id( const_vertex_iterator ) const;
		const_vertex_iterator target_of( const_edge_iterator ) const;
		edge_id_t edge_id( const_edge_iterator ) const;
		bool edge_exists( const_edge_iterator ) const;
		const_edge_iterator edge( const_vertex_iterator, const_vertex_iterator ) const;
		connectivity_address_type edge_value( const_edge_iterator ) const;
		bool edge_exists( const_vertex_iterator, const_vertex_iterator ) const;
		template <std::size_t r_> friend std::ostream & operator<<( std::ostream &, multirelational_graph<r_> const & );
		template <std::size_t r_> friend std::istream & operator>>( std::istream &, multirelational_graph<r_> & );
	private:
		std::size_t num_vertices;
		std::size_t num_edges;
		std::unique_ptr<void*[]> vertices;
		std::unique_ptr<void*[]> edges;
		std::unique_ptr<connectivity_address_type[]> edge_values;
};

template <std::size_t r>
multirelational_graph<r>::multirelational_graph() : num_vertices(0), num_edges(0), vertices(std::unique_ptr<void*[]>(new void*[1])), edges(std::unique_ptr<void*[]>(new void*[1])), edge_values( std::unique_ptr<typename multirelational_graph<r>::connectivity_address_type[]>(new typename multirelational_graph<r>::connectivity_address_type[0])) {
	vertices[0] = &edges[0];
	edges[0] = NULL;
}

template <std::size_t r>
multirelational_graph<r>::multirelational_graph( multirelational_graph const & g ) : num_vertices(g.num_vertices), num_edges(g.num_edges), vertices(std::unique_ptr<void*[]>(new void*[g.vertex_count()+1])), edges(std::unique_ptr<void*[]>(new void*[g.num_edges+1])), edge_values( std::unique_ptr<typename multirelational_graph<r>::connectivity_address_type[]>(new multirelational_graph<r>::connectivity_address_type[g.num_edges])) {
	for( const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it ) {
		vertices[ g.vertex_id( it ) ] = &edges[ g.edge_id( g.neighbors_begin( it ) ) ];
	}
	vertices[ vertex_count() ] = &edges[ num_edges ];
	for( const_edge_iterator it = g.edges_begin(); it != g.edges_end(); ++it ) {
		edges[ g.edge_id( it ) ] = &vertices[ g.vertex_id( g.target_of( it ) ) ];
		edge_values[ g.edge_id( it ) ] = g.edge_value( it );
	}
	edges[ num_edges ] = NULL;
}

template <std::size_t r>
multirelational_graph<r>::~multirelational_graph() {
}

template <std::size_t r>
multirelational_graph<r> & multirelational_graph<r>::operator=( multirelational_graph<r> const & g ) {
	if( this != &g ) {
		num_vertices = g.num_vertices;
		num_edges = g.num_edges;
		vertices = std::unique_ptr<void*[]>(new void*[ g.vertex_count()+1 ]);
		edges = std::unique_ptr<void*[]>(new void*[ g.num_edges+1 ]);
		edge_values = std::unique_ptr<typename multirelational_graph<r>::connectivity_address_type[]>(new multirelational_graph<r>::connectivity_address_type[g.num_edges]);
		for( const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it ) {
			vertex_id_t id = g.vertex_id( it );
			vertices[ id ] = &edges[ g.edge_id( g.neighbors_begin( it ) ) ];
		}
		vertices[ vertex_count() ] = &edges[ num_edges ];
		for( const_edge_iterator it = g.edges_begin(); it != g.edges_end(); ++it ) {
			edge_id_t id = g.edge_id( it );
			edges[ id ] = &vertices[ g.vertex_id( g.target_of( it ) ) ];
			edge_values[ g.edge_id( it ) ] = g.edge_value( it );
		}
		edges[ num_edges ] = NULL;
	}
	return *this;
}

template <std::size_t r>
std::size_t multirelational_graph<r>::vertex_count() const {
	return num_vertices;
}

template <std::size_t r>
std::size_t multirelational_graph<r>::edge_count() const {
	return num_edges / 2;
}

template <std::size_t r>
std::size_t multirelational_graph<r>::relation_count() const {
	static std::size_t count( std::ceil( std::log2( *std::max_element( &edge_values[0], &edge_values[num_edges] ) + 1 ) ) );
	return count;
}

template <std::size_t r>
const_vertex_iterator multirelational_graph<r>::vertices_begin() const {
	return &vertices[ 0 ];
}

template <std::size_t r>
const_vertex_iterator multirelational_graph<r>::vertices_end() const {
	return &vertices[ vertex_count() ];
}

template <std::size_t r>
const_vertex_iterator multirelational_graph<r>::edges_begin() const {
	return &edges[ 0 ];
}

template <std::size_t r>
const_vertex_iterator multirelational_graph<r>::edges_end() const {
	return &edges[ num_edges ];
}

template <std::size_t r>
const_edge_iterator multirelational_graph<r>::neighbors_begin( const_vertex_iterator it ) const {
	return static_cast<const_edge_iterator>( *it );
}

template <std::size_t r>
const_edge_iterator multirelational_graph<r>::neighbors_end( const_vertex_iterator it ) const {
	return static_cast<const_edge_iterator>( *(it+1) );
}

template <std::size_t r>
vertex_id_t multirelational_graph<r>::vertex_id( const_vertex_iterator it ) const {
	return it - static_cast<const_vertex_iterator>( vertices.get() );
}

template <std::size_t r>
const_vertex_iterator multirelational_graph<r>::target_of( const_edge_iterator it ) const {
	return static_cast<const_vertex_iterator>( *it );
}

template <std::size_t r>
edge_id_t multirelational_graph<r>::edge_id( const_edge_iterator it ) const {
	return it - static_cast<const_edge_iterator>( edges.get() );
}

template <std::size_t r>
bool multirelational_graph<r>::edge_exists( const_edge_iterator it ) const {
	return it != edges_end();
}

template <std::size_t r>
const_edge_iterator multirelational_graph<r>::edge( const_vertex_iterator source, const_vertex_iterator target ) const {
	const_edge_iterator it( std::find( neighbors_begin( source ), neighbors_end( source ), target ) );
	return it == neighbors_end( source ) ? edges_end() : it;

}

template <std::size_t r>
typename multirelational_graph<r>::connectivity_address_type multirelational_graph<r>::edge_value( const_edge_iterator it ) const {
	return edge_values[ edge_id( it ) ];
}

template <std::size_t r>
bool multirelational_graph<r>::edge_exists( const_vertex_iterator source, const_vertex_iterator target ) const {
	return neighbors_end( source ) != std::find( neighbors_begin( source ), neighbors_end( source ), target );
}

template <std::size_t r>
std::ostream & operator<<( std::ostream & os, multirelational_graph<r> const & g ) {
	for( const_vertex_iterator vIt = g.vertices_begin(); vIt < g.vertices_end(); ++vIt ) {
		const_edge_iterator nIt = g.neighbors_begin( vIt );
		while( nIt < g.neighbors_end( vIt ) - 1 ) {
			os << g.vertex_id( g.target_of( nIt ) ) << ',' << g.edge_value( nIt ) << ' ';
			++nIt;
		}
		if( nIt < g.neighbors_end( vIt ) ) {
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
std::istream & operator>>( std::istream & is, multirelational_graph<r> & g ) {
	std::vector<edge_id_t> v_temp;
	std::vector<std::pair<vertex_id_t,typename multirelational_graph<r>::connectivity_address_type> > e_temp;
	std::string s1;
	std::size_t edge_id( 0 );
	while( std::getline( is, s1, '\n' ) ) {
		v_temp.push_back( edge_id );
		std::istringstream iss1( s1 );
		std::string s2;
		while( std::getline( iss1, s2, ' ' ) ) {
			std::istringstream iss2( s2 );
			std::string neighbor_str;
			std::string value_str;
			getline( iss2, neighbor_str, ',' );
			getline( iss2, value_str, ',' );
			std::istringstream iss3( neighbor_str );
			std::istringstream iss4( value_str );
			vertex_id_t neighbor;
			typename multirelational_graph<r>::connectivity_address_type value;
			iss3 >> neighbor;
			iss4 >> value;
			e_temp.push_back( std::make_pair( neighbor, value ) );
			++edge_id;
		}
	}
	
	g.num_vertices = v_temp.size();
	g.num_edges = e_temp.size();

	g.vertices = std::unique_ptr<void*[]>(new void*[ g.vertex_count() + 1 ]);
	g.edges = std::unique_ptr<void*[]>(new void*[ g.num_edges + 1 ]);
	g.edge_values = std::unique_ptr<typename multirelational_graph<r>::connectivity_address_type[]>(new typename multirelational_graph<r>::connectivity_address_type[g.num_edges]);

	for( size_t i( 0 ); i < g.vertex_count(); ++i ) {
		g.vertices[ i ] = &g.edges[ v_temp[ i ] ];
	}
	g.vertices[ g.vertex_count() ] = &g.edges[ g.num_edges ];
	for( size_t i( 0 ); i < g.num_edges; ++i ) {
		g.edges[ i ] = &g.vertices[ e_temp[ i ].first ];
		g.edge_values[ i ] = e_temp[ i ].second;
	}
	g.edges[ g.num_edges ] = NULL;
	
	return is;
}

}

#endif
