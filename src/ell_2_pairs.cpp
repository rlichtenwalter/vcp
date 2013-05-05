/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#include <set>
#include <tclap/CmdLine.h>
#include <vcp/graph.hpp>

int main( int argc, char* argv[] ) {
	try {
		TCLAP::CmdLine cmd( "Print, in lexicographical order, all pairs of nodes that are two hops distant from each other.", ' ', "1.0.0" );
		cmd.parse( argc, argv );
	} catch( TCLAP::ArgException & e ) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}
	
	vcp::graph g;
	std::cin >> g;
	for( vcp::const_vertex_iterator vIt( g.vertices_begin() ); vIt != g.vertices_end(); ++vIt ) {
		std::set<vcp::const_vertex_iterator> neighbors;
		for( vcp::const_edge_iterator e1It( g.neighbors_begin( vIt ) ); e1It != g.neighbors_end( vIt ); ++e1It ) {
			for( vcp::const_edge_iterator e2It( g.neighbors_begin( g.target_of( e1It ) ) ); e2It != g.neighbors_end( g.target_of( e1It ) ); ++e2It ) {
				if( vIt < g.target_of( e2It ) && !g.edge_exists( vIt, g.target_of( e2It ) ) ) {
					neighbors.insert( g.target_of( e2It ) );
				}
			}
		}
		for( std::set<vcp::const_vertex_iterator>::const_iterator it( neighbors.begin() ); it != neighbors.end(); ++it ) {
			std::cout << g.vertex_id( vIt ) << ' ' << g.vertex_id( *it ) << '\n';
		}
	}

	return 0;
}
