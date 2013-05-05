/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_VCP_3_R_0
#define VCP_VCP_3_R_0

#include <cstddef>
#include <map>
#include <utility>
#include <vcp/graph.hpp>
#include <vcp/multirelational_graph.hpp>
#include <vcp/vcp_dynamic_mapper.hpp>

namespace vcp {

template <std::size_t n,std::size_t r,bool d> class vcp;

template <std::size_t r>
class vcp<3,r,0> {
	public:
		typedef typename std::conditional< (r>1), multirelational_graph<r>, graph >::type graph_type;
		typedef typename multirelational_graph<r>::connectivity_address_type connectivity_address_type;
		typedef typename vcp_dynamic_mapper<3,r,0>::subgraph_address_type subgraph_address_type;
		vcp( graph_type const & g );
		std::map<subgraph_address_type,unsigned long> const generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 );
	private:
		enum connectivity_value : std::size_t { // in bit shifting terms
			V1V2 = 0 * r,
			V1V3 = 1 * r,
			V2V3 = 2 * r
		};
		graph_type const & g;
};

template <std::size_t r>
vcp<3,r,0>::vcp( graph_type const & g ) : g( g ) {
}

template <std::size_t r>
std::map<typename vcp<3,r,0>::subgraph_address_type,unsigned long> const vcp<3,r,0>::generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 ) {
	std::map<subgraph_address_type,unsigned long> counts;

	subgraph_address_type v1v2( subgraph_address_type( g.edge_value( g.edge( v1, v2 ) ) ) << static_cast<std::size_t>( V1V2 ) );

	const_edge_iterator v1_it( g.neighbors_begin( v1 ) );
	const_edge_iterator v1_end( g.neighbors_end( v1 ) );
	const_edge_iterator v2_it( g.neighbors_begin( v2 ) );
	const_edge_iterator v2_end( g.neighbors_end( v2 ) );
	
	unsigned long union_cardinality( 0 );
	while( v1_it != v1_end && v2_it != v2_end ) {
		if( g.target_of( v1_it ) == v2 ) {
			++v1_it;
		} else if( g.target_of( v2_it ) == v1 ) {
			++v2_it;
		} else {
			++union_cardinality;
			if( g.target_of( v1_it ) < g.target_of( v2_it ) ) {
				++counts.insert( std::make_pair( v1v2 + (subgraph_address_type(g.edge_value( v1_it )) << static_cast<std::size_t>(V1V3)), 0 ) ).first->second;
				++v1_it;
			} else if( g.target_of( v1_it ) > g.target_of( v2_it ) ) {
				++counts.insert( std::make_pair( v1v2 + (subgraph_address_type(g.edge_value( v1_it )) << static_cast<std::size_t>(V2V3)), 0 ) ).first->second;
				++v2_it;
			} else {
				++counts.insert( std::make_pair( v1v2 + (subgraph_address_type(g.edge_value( v1_it )) << static_cast<std::size_t>(V1V3)) + (subgraph_address_type(g.edge_value( v2_it )) <<  static_cast<std::size_t>(V2V3)), 0 ) ).first->second;
				++v1_it;
				++v2_it;
			}
		}
	}
	while( v1_it != v1_end ) {
		if( g.target_of( v1_it ) != v2 ) {
			++union_cardinality;
			++counts.insert( std::make_pair( v1v2 + (subgraph_address_type(g.edge_value( v1_it )) << static_cast<std::size_t>(V1V3)), 0 ) ).first->second;
		}
		++v1_it;
	} while( v2_it != v2_end ) {
		if( g.target_of( v2_it ) != v1 ) {
			++union_cardinality;
			++counts.insert( std::make_pair( v1v2 + (subgraph_address_type(g.edge_value( v2_it )) << static_cast<std::size_t>(V2V3)), 0 ) ).first->second;
		}
		++v2_it;
	}
	counts.insert( std::make_pair( v1v2, g.vertex_count() - 2 - union_cardinality ) );
		
	return counts;
}

}

#endif
