/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_DIRECTED_GRAPH
#define VCP_DIRECTED_GRAPH

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace vcp {

typedef std::size_t vertex_id_t;
typedef std::size_t edge_id_t;
typedef void * const * const_vertex_iterator;
typedef void * const * const_edge_iterator;

class directed_graph {
	public:
		directed_graph();
		directed_graph( directed_graph const & );
		~directed_graph();
		directed_graph & operator=( directed_graph const & );
		std::size_t vertex_count() const;
		std::size_t out_edge_count() const;
		std::size_t in_edge_count() const;
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
		const_edge_iterator out_edge( const_vertex_iterator, const_vertex_iterator ) const;
		const_edge_iterator in_edge( const_vertex_iterator, const_vertex_iterator ) const;
		bool out_edge_exists( const_vertex_iterator, const_vertex_iterator ) const;
		bool in_edge_exists( const_vertex_iterator, const_vertex_iterator ) const;
		friend std::ostream & operator<<( std::ostream &, directed_graph const & );
		friend std::istream & operator>>( std::istream &, directed_graph & );
	private:
		std::size_t num_vertices;
		std::size_t num_out_edges;
		std::unique_ptr<void*[]> vertices;
		std::unique_ptr<void*[]> edges;
};

directed_graph::directed_graph() : num_vertices(0), num_out_edges(0), vertices(std::unique_ptr<void*[]>(new void*[1])), edges(std::unique_ptr<void*[]>(new void*[1])) {
	vertices[0] = &edges[0];
	edges[0] = NULL;
}

directed_graph::directed_graph( directed_graph const & g ) : num_vertices(g.num_vertices), num_out_edges(g.num_out_edges), vertices(std::unique_ptr<void*[]>(new void*[2*g.vertex_count()+1])), edges(std::unique_ptr<void*[]>(new void*[g.out_edge_count()+g.in_edge_count()+1])) {
	for( const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it ) {
		vertices[ g.vertex_id( it ) ] = &edges[ g.edge_id( g.out_neighbors_begin( it ) ) ];
		vertices[ vertex_count() + g.vertex_id( it ) ] = &edges[ g.edge_id( g.in_neighbors_begin( it ) ) ];
	}
	vertices[ 2*vertex_count() ]  = &edges[ out_edge_count() + in_edge_count() ];
	for( const_edge_iterator it = g.out_edges_begin(); it != g.out_edges_end(); ++it ) {
		edges[ g.edge_id( it ) ] = &vertices[ g.vertex_id( g.target_of( it ) ) ];
	}
	for( const_edge_iterator it = g.in_edges_begin(); it != g.in_edges_end(); ++it ) {
		edges[ g.edge_id( it ) ] = &vertices[ g.vertex_id( g.target_of( it ) ) ];
	}
	edges[ out_edge_count() + in_edge_count() ] = NULL;
}

directed_graph::~directed_graph() {
}

directed_graph & directed_graph::operator=( directed_graph const & g ) {
	if( this != &g ) {
		num_vertices = g.num_vertices;
		num_out_edges = g.num_out_edges;
		vertices = std::unique_ptr<void*[]>(new void*[ 2*g.vertex_count()+1 ]);
		edges = std::unique_ptr<void*[]>(new void*[ g.out_edge_count()+g.in_edge_count()+1 ]);
		for( const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it ) {
			vertices[ g.vertex_id( it ) ] = &edges[ g.edge_id( g.out_neighbors_begin( it ) ) ];
			vertices[ vertex_count() + g.vertex_id( it ) ] = &edges[ g.edge_id( g.in_neighbors_begin( it ) ) ];
		}
		vertices[ 2*vertex_count() ]  = &edges[ out_edge_count() + in_edge_count() ];
		for( const_edge_iterator it = g.out_edges_begin(); it != g.out_edges_end(); ++it ) {
			edges[ g.edge_id( it ) ] = &vertices[ g.vertex_id( g.target_of( it ) ) ];
		}
		for( const_edge_iterator it = g.in_edges_begin(); it != g.in_edges_end(); ++it ) {
			edges[ g.edge_id( it ) ] = &vertices[ g.vertex_id( g.target_of( it ) ) ];
		}
		edges[ out_edge_count() + in_edge_count() ] = NULL;
	}
	return *this;
}

std::size_t directed_graph::vertex_count() const {
	return num_vertices;
}

std::size_t directed_graph::out_edge_count() const {
	return num_out_edges;
}

std::size_t directed_graph::in_edge_count() const {
	return num_out_edges;
}

const_vertex_iterator directed_graph::vertices_begin() const {
	return &vertices[ 0 ];
}

const_vertex_iterator directed_graph::vertices_end() const {
	return &vertices[ vertex_count() ];
}

const_vertex_iterator directed_graph::out_edges_begin() const {
	return &edges[ 0 ];
}

const_vertex_iterator directed_graph::out_edges_end() const {
	return &edges[ out_edge_count() ];
}

const_vertex_iterator directed_graph::in_edges_begin() const {
	return &edges[ out_edge_count() ];
}

const_vertex_iterator directed_graph::in_edges_end() const {
	return &edges[ out_edge_count() + in_edge_count() ];
}

const_edge_iterator directed_graph::out_neighbors_begin( const_vertex_iterator it ) const {
	return static_cast<const_edge_iterator>( *it );
}

const_edge_iterator directed_graph::out_neighbors_end( const_vertex_iterator it ) const {
	return static_cast<const_edge_iterator>( *(it+1) );
}

const_edge_iterator directed_graph::in_neighbors_begin( const_vertex_iterator it ) const {
	return static_cast<const_edge_iterator>( *(vertex_count()+it) );
}

const_edge_iterator directed_graph::in_neighbors_end( const_vertex_iterator it ) const {
	return static_cast<const_edge_iterator>( *(vertex_count()+it+1) );
}

vertex_id_t directed_graph::vertex_id( const_vertex_iterator it ) const {
	return it - static_cast<const_vertex_iterator>( vertices.get() );
}

const_vertex_iterator directed_graph::target_of( const_edge_iterator it ) const {
	return static_cast<const_vertex_iterator>( *it );
}

edge_id_t directed_graph::edge_id( const_edge_iterator it ) const {
	return it - static_cast<const_edge_iterator>( edges.get() );
}

bool directed_graph::edge_exists( const_edge_iterator it ) const {
	return it != in_edges_end();
}

const_edge_iterator directed_graph::out_edge( const_vertex_iterator source, const_vertex_iterator target ) const {
	const_edge_iterator it( std::find( out_neighbors_begin( source ), out_neighbors_end( source ), target ) );
	return it == out_neighbors_end( source ) ? in_edges_end() : it;
}

const_edge_iterator directed_graph::in_edge( const_vertex_iterator source, const_vertex_iterator target ) const {
	const_edge_iterator it( std::find( in_neighbors_begin( source ), in_neighbors_end( source ), target ) );
	return it == in_neighbors_end( source ) ? in_edges_end() : it;
}

bool directed_graph::out_edge_exists( const_vertex_iterator source, const_vertex_iterator target ) const {
	return out_neighbors_end( source ) != std::find( out_neighbors_begin( source ), out_neighbors_end( source ), target );
}

bool directed_graph::in_edge_exists( const_vertex_iterator source, const_vertex_iterator target ) const {
	return in_neighbors_end( source ) != std::find( in_neighbors_begin( source ), in_neighbors_end( source ), target );
}

std::ostream & operator<<( std::ostream & os, directed_graph const & g ) {
	for( const_vertex_iterator vIt( g.vertices_begin() ); vIt < g.vertices_end(); ++vIt ) {
		const_edge_iterator nIt( g.out_neighbors_begin( vIt ) );
		while( nIt < g.out_neighbors_end( vIt ) - 1 ) {
			os << g.vertex_id( g.target_of( nIt++ ) ) << ' ';
		}
		if( nIt < g.out_neighbors_end( vIt ) ) {
			os << g.vertex_id( g.target_of( nIt++ ) );
		}
		if( vIt < g.vertices_end() ) {
			os << '\n';
		}
	}
	return os;
}

std::istream & operator>>( std::istream & is, directed_graph & g ) {
	std::vector<edge_id_t> out_v_temp;
	std::vector<vertex_id_t> out_e_temp;
	std::string s1;
	std::size_t edge_id( 0 );
	
	while( std::getline( is, s1, '\n' ) ) {
		out_v_temp.push_back( edge_id );
		std::istringstream iss( s1 );
		std::string s2;
		while( std::getline( iss, s2, ' ' ) ) {
			out_e_temp.push_back( atol( s2.c_str() ) );
			++edge_id;
		}
	}

	g.num_vertices = out_v_temp.size();
	g.num_out_edges = out_e_temp.size();
	
	g.vertices = std::unique_ptr<void*[]>(new void*[ 2 * g.vertex_count() + 1 ]);
	g.edges = std::unique_ptr<void*[]>(new void*[ g.out_edge_count() + g.in_edge_count() + 1 ]);
	
	for( std::size_t i = 0; i < g.vertex_count(); ++i ) {
		g.vertices[ i ] = &g.edges[ out_v_temp.at( i ) ];
	}
	g.vertices[ g.vertex_count() ] = &g.edges[ g.out_edge_count() ];
	for( std::size_t i = 0; i < g.out_edge_count(); ++i ) {
		g.edges[ i ] = &g.vertices[ out_e_temp.at( i ) ];
	}

	out_v_temp = std::vector<edge_id_t>();
	out_e_temp = std::vector<vertex_id_t>();

	std::vector<std::vector<vertex_id_t> > in_temp( g.vertex_count() );
	for( const_vertex_iterator vIt( g.vertices_begin() ); vIt != g.vertices_end(); ++vIt ) {
		for( const_edge_iterator eIt( g.out_neighbors_begin( vIt ) ); eIt != g.out_neighbors_end( vIt ); ++eIt ) {
			in_temp[ g.vertex_id( g.target_of( eIt ) ) ].push_back( g.vertex_id( vIt ) );
		}
	}

	std::size_t in_count( 0 );
	for( std::size_t i = 0; i < g.vertex_count(); ++i ) {
		g.vertices[ g.vertex_count()+i ] = &g.edges[ g.out_edge_count() + in_count ];
		std::vector<vertex_id_t> const & in_neighbors = in_temp[ i ];
		for( std::vector<vertex_id_t>::const_iterator it( in_neighbors.begin() ); it != in_neighbors.end(); ++it ) {
			g.edges[ g.out_edge_count() + in_count++ ] = &g.vertices[ *it ];
		}
	}

	g.vertices[ 2 * g.vertex_count() ] = &g.edges[ g.out_edge_count() + g.in_edge_count() ];
	g.edges[ g.out_edge_count() + g.in_edge_count() ] = NULL;
	
	return is;
}

}

#endif
