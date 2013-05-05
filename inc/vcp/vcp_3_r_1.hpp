/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_VCP_3_R_1
#define VCP_VCP_3_R_1

#include <cstddef>
#include <map>
#include <utility>
#include <vcp/directed_graph.hpp>
#include <vcp/multirelational_directed_graph.hpp>
#include <vcp/vcp_dynamic_mapper.hpp>

namespace vcp {
	
template <std::size_t n,std::size_t r,bool d> class vcp;

template <std::size_t r>
class vcp<3,r,1> {
	public:
		typedef typename std::conditional< (r>1), multirelational_directed_graph<r>, directed_graph >::type graph_type;
		typedef typename multirelational_graph<r>::connectivity_address_type connectivity_address_type;
		typedef typename vcp_dynamic_mapper<3,r,1>::subgraph_address_type subgraph_address_type;
		vcp( graph_type const & g );
		std::map<subgraph_address_type,unsigned long> const generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 );
	private:
		enum connectivity_value : std::size_t { // in bit shifting terms
			V1V2 = 0 * r,
			V1V3 = 2 * r,
			V2V3 = 4 * r
		};
		enum directedness_value : std::size_t { // in bit shifting terms
			OUT = 0 * r,
			IN = 1 * r
		};
		graph_type const & g;
		std::pair<const_edge_iterator,std::pair<connectivity_address_type,connectivity_address_type> > next_union_element( const_edge_iterator & it1, const_edge_iterator end1, const_edge_iterator & it2, const_edge_iterator end2 ) const;
};

template <std::size_t r>
vcp<3,r,1>::vcp( graph_type const & g ) : g(g) {
}

template <std::size_t r>
std::pair<const_edge_iterator,std::pair<typename vcp<3,r,1>::connectivity_address_type,typename vcp<3,r,1>::connectivity_address_type> > vcp<3,r,1>::next_union_element( const_edge_iterator & it1, const_edge_iterator end1, const_edge_iterator & it2, const_edge_iterator end2 ) const { // out-neighbor iterators should always come first
	std::pair<const_edge_iterator,std::pair<connectivity_address_type,connectivity_address_type> > temp;
	if( it1 != end1 && it2 != end2 ) {
		if( g.target_of( it1 ) < g.target_of( it2 ) ) {
			temp = std::make_pair( it1, std::make_pair( g.edge_value( it1 ), 0 ) );
			++it1;
		} else if( g.target_of( it1 ) > g.target_of( it2 ) ) {
			temp = std::make_pair( it2, std::make_pair( 0, g.edge_value( it2 ) ) );
			++it2;
		} else {
			temp = std::make_pair( it1, std::make_pair( g.edge_value( it1 ), g.edge_value( it2 ) ) );
			++it1;
			++it2;
		}
	} else if( it1 != end1 ) {
			temp = std::make_pair( it1, std::make_pair( g.edge_value( it1 ), 0 ) );
			++it1;
	} else if( it2 != end2 ) {
			temp = std::make_pair( it2, std::make_pair( 0, g.edge_value( it2 ) ) );
			++it2;
	} else {
		return std::make_pair( end2, std::make_pair( 0, 0 ) ); // the second element of this pair should never be used
	}
	return temp;
}

template <std::size_t r>
std::map<typename vcp<3,r,1>::subgraph_address_type,unsigned long> const vcp<3,r,1>::generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 ) {
	std::map<subgraph_address_type,unsigned long> counts;
	
	subgraph_address_type v1v2( (subgraph_address_type( g.edge_value( g.out_edge( v1, v2 ) ) ) << (V1V2 + OUT)) + (subgraph_address_type( g.edge_value( g.in_edge( v1, v2 ) ) ) << (V1V2 + IN )) );
	
	const_edge_iterator v1_out_neighbors_it( g.out_neighbors_begin( v1 ) );
	const_edge_iterator v1_out_neighbors_end( g.out_neighbors_end( v1 ) );
	const_edge_iterator v1_in_neighbors_it( g.in_neighbors_begin( v1 ) );
	const_edge_iterator v1_in_neighbors_end( g.in_neighbors_end( v1 ) );
	const_edge_iterator v2_out_neighbors_it( g.out_neighbors_begin( v2 ) );
	const_edge_iterator v2_out_neighbors_end( g.out_neighbors_end( v2 ) );
	const_edge_iterator v2_in_neighbors_it( g.in_neighbors_begin( v2 ) );
	const_edge_iterator v2_in_neighbors_end( g.in_neighbors_end( v2 ) );
	std::pair<const_edge_iterator,std::pair<connectivity_address_type,connectivity_address_type> > min1( next_union_element( v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it, v1_in_neighbors_end ) );
	std::pair<const_edge_iterator,std::pair<connectivity_address_type,connectivity_address_type> > min2( next_union_element( v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it, v2_in_neighbors_end ) );
	unsigned long union_cardinality( 0 );
	while( min1.first != v1_in_neighbors_end && min2.first != v2_in_neighbors_end ) {
		if( g.target_of( min1.first ) < g.target_of( min2.first ) ) {
			if( g.target_of( min1.first ) != v2 ) {
				++union_cardinality;
				++counts.insert( std::make_pair( v1v2 + (subgraph_address_type( min1.second.first ) << (V1V3 + OUT)) + (subgraph_address_type( min1.second.second ) << (V1V3 + IN)), 0 ) ).first->second;
			}
			min1 = next_union_element( v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it, v1_in_neighbors_end );
		} else if( g.target_of( min1.first ) > g.target_of( min2.first ) ) {
			if( g.target_of( min2.first ) != v1 ) {
				++union_cardinality;
				++counts.insert( std::make_pair( v1v2 + (subgraph_address_type( min2.second.first ) << (V2V3 + OUT)) + (subgraph_address_type( min2.second.second ) << (V2V3 + IN)), 0 ) ).first->second;
			}
			min2 = next_union_element( v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it, v2_in_neighbors_end );
		} else { // the next neighbor is shared by both v1 and v2, so it cannot be either and we do not need to check to exclude it
			++union_cardinality;
			++counts.insert( std::make_pair( v1v2 + ((subgraph_address_type( min1.second.first )) << (V1V3 + OUT)) + (subgraph_address_type( min1.second.second ) << (V1V3 + IN)) + (subgraph_address_type( min2.second.first ) << (V2V3 + OUT)) + (subgraph_address_type( min2.second.second ) << (V2V3 + IN)), 0 ) ).first->second;
			min1 = next_union_element( v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it, v1_in_neighbors_end );
			min2 = next_union_element( v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it, v2_in_neighbors_end );
		}
	}
	while( min1.first != v1_in_neighbors_end ) {
		if( g.target_of( min1.first ) != v2 ) {
			++union_cardinality;
			++counts.insert( std::make_pair( v1v2 + (subgraph_address_type( min1.second.first ) << (V1V3 + OUT)) + (subgraph_address_type( min1.second.second ) << (V1V3 + IN)), 0 ) ).first->second;
		}
		min1 = next_union_element( v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it, v1_in_neighbors_end );
	}
	while( min2.first != v2_in_neighbors_end ) {
		if( g.target_of( min2.first ) != v1 ) {
			++union_cardinality;
			++counts.insert( std::make_pair( v1v2 + (subgraph_address_type( min2.second.first ) << (V2V3 + OUT)) + (subgraph_address_type( min2.second.second ) << (V2V3 + IN)), 0 ) ).first->second;
		}
		min2 = next_union_element( v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it, v2_in_neighbors_end );
	}
	
	counts.insert( std::make_pair( v1v2, g.vertex_count() - 2 - union_cardinality ) );

	return counts;
}

}

#endif
