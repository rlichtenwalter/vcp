/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_GRAPH_H
#define VCP_GRAPH_H

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

class graph {
	public:
		graph();
		graph( graph const & );
		graph & operator=( graph const & );
		~graph();
		std::size_t vertex_count() const;
		std::size_t edge_count() const;
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
		bool edge_exists( const_vertex_iterator, const_vertex_iterator ) const;
		friend std::ostream & operator<<( std::ostream &, graph const & );
		friend std::istream & operator>>( std::istream &, graph & );
	private:
		std::size_t num_vertices;
		std::size_t num_edges;
		std::unique_ptr<void*[]> vertices;
		std::unique_ptr<void*[]> edges;
};

graph::graph() : num_vertices(0), num_edges(0), vertices(std::unique_ptr<void*[]>(new void*[1])), edges(std::unique_ptr<void*[]>(new void*[1])) {
	vertices[0] = &edges[0];
	edges[0] = NULL;
}

graph::graph( graph const & g ) : num_vertices(g.num_vertices), num_edges(g.num_edges), vertices(std::unique_ptr<void*[]>(new void*[g.vertex_count()+1])), edges(std::unique_ptr<void*[]>(new void*[g.num_edges+1])) {
	for( const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it ) {
		vertices[ g.vertex_id( it ) ] = &edges[ g.edge_id( g.neighbors_begin( it ) ) ];
	}
	vertices[ vertex_count() ] = &edges[ num_edges ];
	for( const_edge_iterator it = g.edges_begin(); it != g.edges_end(); ++it ) {
		edges[ g.edge_id( it ) ] = &vertices[ g.vertex_id( g.target_of( it ) ) ];
	}
	edges[ num_edges ] = NULL;
}

graph::~graph() {
}

graph & graph::operator=( graph const & g ) {
	if( this != &g ) {
		num_vertices = g.num_vertices;
		num_edges = g.num_edges;
		vertices = std::unique_ptr<void*[]>(new void*[ g.vertex_count()+1 ]);
		edges = std::unique_ptr<void*[]>(new void*[ g.num_edges+1 ]);
		for( const_vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it ) {
			vertex_id_t id = g.vertex_id( it );
			vertices[ id ] = &edges[ g.edge_id( g.neighbors_begin( it ) ) ];
		}
		vertices[ vertex_count() ] = &edges[ num_edges ];
		for( const_edge_iterator it = g.edges_begin(); it != g.edges_end(); ++it ) {
			edge_id_t id = g.edge_id( it );
			edges[ id ] = &vertices[ g.vertex_id( g.target_of( it ) ) ];
		}
		edges[ num_edges ] = NULL;
	}
	return *this;
}

std::size_t graph::vertex_count() const {
	return num_vertices;
}

std::size_t graph::edge_count() const {
	return num_edges / 2;
}

const_vertex_iterator graph::vertices_begin() const {
	return &vertices[ 0 ];
}

const_vertex_iterator graph::vertices_end() const {
	return &vertices[ vertex_count() ];
}

const_vertex_iterator graph::edges_begin() const {
	return &edges[ 0 ];
}

const_vertex_iterator graph::edges_end() const {
	return &edges[ num_edges ];
}

const_edge_iterator graph::neighbors_begin( const_vertex_iterator it ) const {
	return static_cast<const_edge_iterator>( *it );
}

const_edge_iterator graph::neighbors_end( const_vertex_iterator it ) const {
	return static_cast<const_edge_iterator>( *(it+1) );
}

vertex_id_t graph::vertex_id( const_vertex_iterator it ) const {
	return it - static_cast<const_vertex_iterator>( vertices.get() );
}

const_vertex_iterator graph::target_of( const_edge_iterator it ) const {
	return static_cast<const_vertex_iterator>( *it );
}

edge_id_t graph::edge_id( const_edge_iterator it ) const {
	return it - static_cast<const_edge_iterator>( edges.get() );
}

bool graph::edge_exists( const_edge_iterator it ) const {
	return it != edges_end();
}

const_edge_iterator graph::edge( const_vertex_iterator source, const_vertex_iterator target ) const {
	const_edge_iterator it( std::find( neighbors_begin( source ), neighbors_end( source ), target ) );
	return it == neighbors_end( source ) ? edges_end() : it;
}

bool graph::edge_exists( const_vertex_iterator source, const_vertex_iterator target ) const {
	return neighbors_end( source ) != std::find( neighbors_begin( source ), neighbors_end( source ), target );
}

std::ostream & operator<<( std::ostream & os, graph const & g ) {
	for( const_vertex_iterator vIt = g.vertices_begin(); vIt < g.vertices_end(); ++vIt ) {
		const_edge_iterator nIt = g.neighbors_begin( vIt );
		while( nIt < g.neighbors_end( vIt ) - 1 ) {
			os << g.vertex_id( g.target_of( nIt++ ) ) << ' ';
		}
		if( nIt < g.neighbors_end( vIt ) ) {
			os << g.vertex_id( g.target_of( nIt++ ) );
		}
		if( vIt < g.vertices_end() ) {
			os << '\n';
		}
	}
	return os;
}

std::istream & operator>>( std::istream & is, graph & g ) {
	std::vector<edge_id_t> v_temp;
	std::vector<vertex_id_t> e_temp;
	std::string s1;
	std::size_t edge_id( 0 );
	while( std::getline( is, s1, '\n' ) ) {
		v_temp.push_back( edge_id );
		std::istringstream iss( s1 );
		std::string s2;
		while( std::getline( iss, s2, ' ' ) ) {
			e_temp.push_back( atol( s2.c_str() ) );
			++edge_id;
		}
	}
	
	g.num_vertices = v_temp.size();
	g.num_edges = e_temp.size();

	g.vertices = std::unique_ptr<void*[]>(new void*[ g.vertex_count() + 1 ]);
	g.edges = std::unique_ptr<void*[]>(new void*[ g.num_edges + 1 ]);

	for( size_t i = 0; i < g.vertex_count(); ++i ) {
		g.vertices[ i ] = &g.edges[ v_temp.at( i ) ];
	}
	g.vertices[ g.vertex_count() ] = &g.edges[ g.num_edges ];
	for( size_t i = 0; i < g.num_edges; ++i ) {
		g.edges[ i ] = &g.vertices[ e_temp.at( i ) ];
	}
	g.edges[ g.num_edges ] = NULL;
	
	return is;
}

}

#endif
