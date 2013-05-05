/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>
#include <tclap/CmdLine.h>
#include <vcp/directed_graph.hpp>

int main( int argc, char* argv[] ) {
	bool bidirectional;
	try {
		TCLAP::CmdLine cmd( "Convert a directed graph to an undirected graph by considering all unidirectional edges bidirectional.", ' ', "1.0.0" );
		TCLAP::SwitchArg bidirectionalArg( "b", "bidirectional", "Only bidirectional edges become undirected edges in the output. Unidirectional edges are removed.", cmd );
		cmd.parse( argc, argv );
		bidirectional = bidirectionalArg.isSet();
	} catch( TCLAP::ArgException & e ) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}
	
	vcp::directed_graph g;
	std::cin >> g;
	std::vector<vcp::const_vertex_iterator> neighbors;
	for( vcp::const_vertex_iterator vIt( g.vertices_begin() ); vIt != g.vertices_end(); ++vIt ) {
		vcp::const_edge_iterator outIt = g.out_neighbors_begin( vIt );
		vcp::const_edge_iterator outEnd = g.out_neighbors_end( vIt );
		vcp::const_edge_iterator inIt = g.in_neighbors_begin( vIt );
		vcp::const_edge_iterator inEnd = g.in_neighbors_end( vIt );
		while( outIt != outEnd && inIt != inEnd ) {
			if( g.target_of( outIt ) < g.target_of( inIt ) ) {
				if( !bidirectional ) {
					neighbors.push_back( g.target_of( outIt ) );
				}
				++outIt;
			} else if( g.target_of( outIt ) > g.target_of( inIt ) ) {
				if( !bidirectional ) {
					neighbors.push_back( g.target_of( inIt ) );
				}
				++inIt;
			} else {
				neighbors.push_back( g.target_of( outIt ) );
				++outIt;
				++inIt;
			}
		}
		if( !bidirectional ) {
			while( outIt != outEnd ) {
				neighbors.push_back( g.target_of( outIt ) );
			}
			while( inIt != inEnd ) {
				neighbors.push_back( g.target_of( inIt ) );
			}
		}
		std::vector<vcp::const_vertex_iterator>::const_iterator neighbors_it( neighbors.begin() );
		while( neighbors_it < neighbors.end() - 1 ) {
			std::cout << g.vertex_id( *neighbors_it++ ) << ' ';
		}
		if( neighbors_it < neighbors.end() ) {
			std::cout << g.vertex_id( *neighbors_it );
		}
		std::cout << '\n';
		neighbors.clear();
	}

	return 0;
}
