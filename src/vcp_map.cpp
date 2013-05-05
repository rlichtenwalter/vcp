/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <tclap/CmdLine.h>
#include <vcp/vcp_static_mapper.hpp>

class MemoryArgConstraint : public TCLAP::Constraint<std::string> {
	public:
		std::string description() const {
			return "String must be entirely numeric or numeric ended by [kKmMgG] (e.g. 1024,500m,2G,4g).";
		}
		std::string shortID() const {
			return "MAXIMUM_MEMORY";
		}
		bool check( std::string const & value ) const {
			if( value.empty() ) {
				return false;
			}
			std::string::const_iterator it( value.begin() );
			while( it != value.end() - 1 ) {
				if( !std::isdigit( *it ) ) {
					return false;
				}
				++it;
			}
			if( !std::isdigit( *it ) ) {
				if( value.size() < 2 || (std::tolower( *it ) != 'k' && std::tolower( *it ) != 'm' && std::tolower( *it ) != 'g') ) {
					return false;
				}
			}
			return true;
		}
};

int main( int argc, const char* argv[] ) {
	std::size_t n;
	std::size_t r;
	bool d;
	std::size_t max_bytes( std::numeric_limits<std::size_t>::max() );
	try {
		TCLAP::CmdLine cmd( "Output the subgraph-to-element mapping for a particular VCP", ' ', "1.0.0" );
		MemoryArgConstraint mac;
		TCLAP::ValueArg<std::string> memoryArg( "m", "mmax", "If more than MAXIMUM_MEMORY would be required for computation, do not proceed with computation. Print the number of subgraphs in the VCP and exit.", false, "", &mac, cmd );
		std::vector<std::size_t> allowedN {3, 4, 5, 6, 7, 8};
		TCLAP::ValuesConstraint<std::size_t> allowedNVals( allowedN );
		TCLAP::UnlabeledValueArg<std::size_t> nArg( "n", "n\tNumber of vertices in the VCP", true, 3, &allowedNVals, cmd );
		std::vector<std::size_t> allowedR {1, 2, 3, 4, 5, 6, 7, 8};
		TCLAP::ValuesConstraint<std::size_t> allowedRVals( allowedR );
		TCLAP::UnlabeledValueArg<std::size_t> rArg( "r", "r\tNumber of relations in the VCP", true, 1, &allowedRVals, cmd );
		std::vector<std::size_t> allowedD {0, 1};
		TCLAP::ValuesConstraint<std::size_t> allowedDVals( allowedD );
		TCLAP::UnlabeledValueArg<std::size_t> dArg( "d", "d\tWhether the VCP considers directedness", true, 0, &allowedDVals, cmd );
		cmd.parse( argc, argv );
		n = nArg.getValue();
		r = rArg.getValue();
		d = dArg.getValue();
		std::string s( memoryArg.getValue() );
		if( !s.empty() ) {
			max_bytes = 1;
			char suffix( std::tolower( *(s.end() - 1) ) );
			if( suffix == 'k' ) {
				s = s.substr( 0, s.size() - 1 );
				max_bytes *= 1024;
			} else if( suffix == 'm' ) {
				s = s.substr( 0, s.size() - 11 );
				max_bytes *= 1024 * 1024;
			} else if( suffix == 'g' ) {
				s = s.substr( 0, s.size() - 1 );
				max_bytes *= 1024 * 1024 * 1024;
			}
			std::size_t val;
			std::stringstream( s ) >> val;
			max_bytes *= val;
		}
	} catch( TCLAP::ArgException & e ) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}
	
	std::size_t subgraph_count( vcp::vcp_static_mapper::subgraph_count( n, r, d ) );
	if( subgraph_count * sizeof( std::size_t ) > max_bytes ) {
		std::cout << subgraph_count << std::endl;
		return 0;
	}
	vcp::vcp_static_mapper mapper( n, r, d );
	std::cout << mapper;
	
	return 0;
}

