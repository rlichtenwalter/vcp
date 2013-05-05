/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <tclap/CmdLine.h>
#include <vcp/directed_graph.hpp>
#include <vcp/graph.hpp>
#include <vcp/multirelational_graph.hpp>
#include <vcp/multirelational_directed_graph.hpp>
#include <vcp/vcp.hpp>

template <std::size_t n>
std::ostream & operator<<( std::ostream & os, std::array<unsigned long,n> const & array ) {
	if( array.empty() ) {
		return os;
	}
	std::cout << array[0];
	for( std::size_t i( 1 ); i < n; ++i ) {
		std::cout << ' ' << array[i];
	}
	std::cout << '\n';
	return os;
}

template <typename T>
std::ostream & operator<<( std::ostream & os, std::map<T,unsigned long> const & map ) {
	if( map.empty() ) {
		return os;
	}
	typename std::map<T,unsigned long>::const_iterator it( map.begin() );
	std::cout << it->first << ',' << it->second;
	++it;
	for( ; it != map.end(); ++it ) {
		std::cout << ' ' << it->first << ',' << it->second;
	}
	std::cout << '\n';
	return os;
}

int main( int argc, char * argv[] ) {
	std::size_t n;
	std::size_t r;
	bool d;
	std::string filename;
	try {
		TCLAP::CmdLine cmd( "Output VCP vectors for pairs read from standard input.", ' ', "1.0.0" );
		std::vector<std::size_t> allowedN {3, 4, 5, 6, 7, 8};
		TCLAP::ValuesConstraint<std::size_t> allowedNVals( allowedN );
		TCLAP::UnlabeledValueArg<std::size_t> nArg( "n", "n\tNumber of vertices in the VCP", true, 3, &allowedNVals, cmd );
		TCLAP::UnlabeledValueArg<std::size_t> rArg( "r", "r\tNumber of relations in the VCP", true, 1, "[1,inf]", cmd );
		std::vector<std::size_t> allowedD {0, 1};
		TCLAP::ValuesConstraint<std::size_t> allowedDVals( allowedD );
		TCLAP::UnlabeledValueArg<std::size_t> dArg( "d", "d\tWhether the VCP considers directedness", true, 0, &allowedDVals, cmd );
		TCLAP::UnlabeledValueArg<std::string> filenameArg( "graph_filename", "\tThe name of the file containing the graph", true, "", "graph_filename", cmd );
		cmd.parse( argc, argv );
		n = nArg.getValue();
		r = rArg.getValue();
		d = dArg.getValue();
		filename = filenameArg.getValue();
	} catch( TCLAP::ArgException & e ) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}
	
	std::ifstream file;
	file.open( filename, std::ifstream::in );
	if( !file ) {
		std::cerr << "error opening file: " << filename << std::endl;
	}

	vcp::vertex_id_t v1;
	vcp::vertex_id_t v2;

	if( d ) {
		if( n == 3 ) {
			if( r == 1 ) {
				vcp::directed_graph g;
				file >> g;
				vcp::vcp<3,1,1> profiler( g );
				while( std::cin >> v1 >> v2 ) {
					std::cout << profiler.generate_vector( vcp::const_vertex_iterator( g.vertices_begin() + v1 ), vcp::const_vertex_iterator( g.vertices_begin() + v2) );
				}
			} else if( r == 2 ) {
				vcp::multirelational_directed_graph<2> g;
				file >> g;
				vcp::vcp<3,2,1> profiler( g );
				while( std::cin >> v1 >> v2 ) {
					std::cout << profiler.generate_vector( vcp::const_vertex_iterator( g.vertices_begin() + v1 ), vcp::const_vertex_iterator( g.vertices_begin() + v2) );
				}
			} else if( r == 30 ) {
				vcp::multirelational_directed_graph<30> g;
				file >> g;
				vcp::vcp<3,30,1> profiler( g );
				while( std::cin >> v1 >> v2 ) {
					std::cout << profiler.generate_vector( vcp::const_vertex_iterator( g.vertices_begin() + v1 ), vcp::const_vertex_iterator( g.vertices_begin() + v2) );
				}
			}
		} else if( n == 4 ) {
			if( r == 1 ) {
				vcp::directed_graph g;
				file >> g;
				vcp::vcp<4,1,1> profiler( g );
				while( std::cin >> v1 >> v2 ) {
					std::cout << profiler.generate_vector( vcp::const_vertex_iterator( g.vertices_begin() + v1 ), vcp::const_vertex_iterator( g.vertices_begin() + v2) );
				}
			} else if( r == 2 ) {
				vcp::multirelational_directed_graph<2> g;
				file >> g;
				vcp::vcp<4,2,1> profiler( g );
				while( std::cin >> v1 >> v2 ) {
					std::cout << profiler.generate_vector( vcp::const_vertex_iterator( g.vertices_begin() + v1 ), vcp::const_vertex_iterator( g.vertices_begin() + v2) );
				}
			} else if( r == 30 ) {
				vcp::multirelational_directed_graph<30> g;
				file >> g;
				vcp::vcp<4,30,1> profiler( g );
				while( std::cin >> v1 >> v2 ) {
					std::cout << profiler.generate_vector( vcp::const_vertex_iterator( g.vertices_begin() + v1 ), vcp::const_vertex_iterator( g.vertices_begin() + v2) );
				}
			}
		}
	} else {
		if( n == 3 ) {
			if( r == 1 ) {
				vcp::graph g;
				file >> g;
				vcp::vcp<3,1,0> profiler( g );
				while( std::cin >> v1 >> v2 ) {
					std::cout << profiler.generate_vector( vcp::const_vertex_iterator( g.vertices_begin() + v1 ), vcp::const_vertex_iterator( g.vertices_begin() + v2) );
				}
			} else if( r == 2 ) {
				vcp::multirelational_graph<2> g;
				file >> g;
				vcp::vcp<3,2,0> profiler( g );
				while( std::cin >> v1 >> v2 ) {
					std::cout << profiler.generate_vector( vcp::const_vertex_iterator( g.vertices_begin() + v1 ), vcp::const_vertex_iterator( g.vertices_begin() + v2) );
				}
			} else if( r == 30 ) {
				vcp::multirelational_graph<30> g;
				file >> g;
				vcp::vcp<3,30,0> profiler( g );
				while( std::cin >> v1 >> v2 ) {
					std::cout << profiler.generate_vector( vcp::const_vertex_iterator( g.vertices_begin() + v1 ), vcp::const_vertex_iterator( g.vertices_begin() + v2) );
				}
			}
		} else if( n == 4 ) {
			if( r == 1 ) {
				vcp::graph g;
				file >> g;
				vcp::vcp<4,1,0> profiler( g );
				while( std::cin >> v1 >> v2 ) {
					std::cout << profiler.generate_vector( vcp::const_vertex_iterator( g.vertices_begin() + v1 ), vcp::const_vertex_iterator( g.vertices_begin() + v2) );
				}
			} else if( r == 2 ) {
				vcp::multirelational_graph<2> g;
				file >> g;
				vcp::vcp<4,2,0> profiler( g );
				while( std::cin >> v1 >> v2 ) {
					std::cout << profiler.generate_vector( vcp::const_vertex_iterator( g.vertices_begin() + v1 ), vcp::const_vertex_iterator( g.vertices_begin() + v2) );
				}
			} else if( r == 30 ) {
				vcp::multirelational_graph<30> g;
				file >> g;
				vcp::vcp<4,30,0> profiler( g );
				while( std::cin >> v1 >> v2 ) {
					std::cout << profiler.generate_vector( vcp::const_vertex_iterator( g.vertices_begin() + v1 ), vcp::const_vertex_iterator( g.vertices_begin() + v2) );
				}
			}
		}
	}

	return 0;
}
