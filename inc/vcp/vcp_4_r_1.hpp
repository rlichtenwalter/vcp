/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_VCP_4_R_1
#define VCP_VCP_4_R_1

#include <cassert>
#include <cstddef>
#include <map>
#include <utility>
#include <vcp/multirelational_directed_graph.hpp>
#include <vcp/square_matrix.hpp>
#include <vcp/vcp_dynamic_mapper.hpp>

namespace vcp {

template <std::size_t n,std::size_t r,bool d> class vcp;

template <std::size_t r>
class vcp<4,r,1> {
	public:
		typedef typename multirelational_directed_graph<r>::connectivity_address_type connectivity_address_type;
		typedef typename vcp_dynamic_mapper<4,r,1>::subgraph_address_type subgraph_address_type;
		vcp( multirelational_directed_graph<r> const & );
		std::map<subgraph_address_type,unsigned long> const generate_vector( const_vertex_iterator, const_vertex_iterator );
	private:
		typedef square_matrix<connectivity_address_type,4> connectivity_matrix;
		multirelational_directed_graph<r> const & g;
		vcp_dynamic_mapper<4,r,1> mapper;
		std::map<std::pair<connectivity_address_type,connectivity_address_type>,unsigned long> edge_types;
		std::unique_ptr<std::pair<const_vertex_iterator,connectivity_matrix >[]> v3Vertices;
		std::pair<const_edge_iterator,std::pair<connectivity_address_type,connectivity_address_type> > next_union_element( const_edge_iterator &, const_edge_iterator, const_edge_iterator &, const_edge_iterator ) const;

};

template <std::size_t r>
vcp<4,r,1>::vcp( multirelational_directed_graph<r> const & g ) : g(g), mapper(), v3Vertices( std::unique_ptr<std::pair<const_vertex_iterator,connectivity_matrix>[]>(new std::pair<const_vertex_iterator,connectivity_matrix>[ MAX_NEIGHBORS ] )) {
	unsigned long & gaps( edge_types.insert( std::make_pair( std::make_pair( 0, 0 ), g.vertex_count() * (g.vertex_count() - 1) / 2 ) ).first->second );
	for( const_vertex_iterator it( g.vertices_begin() ); it != g.vertices_end(); ++it ) {
		const_edge_iterator outIt = g.out_neighbors_begin( it );
		const_edge_iterator outEnd = g.out_neighbors_end( it );
		const_edge_iterator inIt = g.in_neighbors_begin( it );
		const_edge_iterator inEnd = g.in_neighbors_end( it );
		while( outIt != outEnd && g.target_of( outIt ) <= it ) {
			++outIt;
		}
		while( inIt != inEnd && g.target_of( inIt ) <= it ) {
			++inIt;
		}
		while( outIt != outEnd && inIt != inEnd ) {
			if( g.target_of( outIt ) < g.target_of( inIt ) ) {
				++edge_types.insert( std::make_pair( std::make_pair( 0, g.edge_value( outIt ) ), 0 ) ).first->second;
				--gaps;
				++outIt;
			} else if( g.target_of( outIt ) > g.target_of( inIt ) ) {
				++edge_types.insert( std::make_pair( std::make_pair( 0, g.edge_value( inIt ) ), 0 ) ).first->second;
				--gaps;
				++inIt;
			} else {
				++edge_types.insert( std::make_pair( g.edge_value( outIt ) < g.edge_value( inIt ) ? std::make_pair( g.edge_value( outIt ), g.edge_value( inIt ) ) : std::make_pair( g.edge_value( inIt ), g.edge_value( outIt ) ), 0 ) ).first->second;
				--gaps;
				++outIt;
				++inIt;
			}
		}
		while( outIt != outEnd ) {
			++edge_types.insert( std::make_pair( std::make_pair( 0, g.edge_value( outIt ) ), 0 ) ).first->second;
			--gaps;
			++outIt;
		}
		while( inIt != inEnd ) {
			++edge_types.insert( std::make_pair( std::make_pair( 0, g.edge_value( inIt ) ), 0 ) ).first->second;
			--gaps;
			++inIt;
		}		
	}
}

template <std::size_t r>
std::pair<const_edge_iterator,std::pair<typename vcp<4,r,1>::connectivity_address_type,typename vcp<4,r,1>::connectivity_address_type> > vcp<4,r,1>::next_union_element( const_edge_iterator & it1, const_edge_iterator end1, const_edge_iterator & it2, const_edge_iterator end2 ) const {
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
std::map<typename vcp<4,r,1>::subgraph_address_type,unsigned long> const vcp<4,r,1>::generate_vector( const_vertex_iterator v1, const_vertex_iterator v2 ) {
	std::map<subgraph_address_type,unsigned long> counts;
	std::map<std::pair<connectivity_address_type,connectivity_address_type>,unsigned long> temp_edge_types;
	
	connectivity_matrix connectivity;
	connectivity( 0, 1 ) = g.edge_value( g.out_edge( v1, v2 ) );
	connectivity( 1, 0 ) = g.edge_value( g.in_edge( v1, v2 ) );
	
	unsigned long & gaps( temp_edge_types.insert( std::make_pair( std::make_pair( 0, 0 ), 0 ) ).first->second );
	
	// compose ordered list of v3 candidates
	const_edge_iterator v1_out_neighbors_it( g.out_neighbors_begin( v1 ) );
	const_edge_iterator v1_out_neighbors_end( g.out_neighbors_end( v1 ) );
	const_edge_iterator v1_in_neighbors_it( g.in_neighbors_begin( v1 ) );
	const_edge_iterator v1_in_neighbors_end( g.in_neighbors_end( v1 ) );
	const_edge_iterator v2_out_neighbors_it( g.out_neighbors_begin( v2 ) );
	const_edge_iterator v2_out_neighbors_end( g.out_neighbors_end( v2 ) );
	const_edge_iterator v2_in_neighbors_it( g.in_neighbors_begin( v2 ) );
	const_edge_iterator v2_in_neighbors_end( g.in_neighbors_end( v2 ) );
	assert( MAX_NEIGHBORS > (v1_out_neighbors_end-v1_out_neighbors_it)+(v1_in_neighbors_end-v1_in_neighbors_it)+(v2_out_neighbors_end-v2_out_neighbors_it)+(v2_in_neighbors_end-v2_in_neighbors_it) );
	std::pair<const_vertex_iterator,connectivity_matrix>* v3Vertices_begin( &v3Vertices[0] );
	std::pair<const_vertex_iterator,connectivity_matrix>* v3Vertices_end( &v3Vertices[0] );
	std::pair<const_edge_iterator,std::pair<connectivity_address_type,connectivity_address_type> > min1( next_union_element( v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it, v1_in_neighbors_end ) );
	std::pair<const_edge_iterator,std::pair<connectivity_address_type,connectivity_address_type> > min2( next_union_element( v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it, v2_in_neighbors_end ) );
	while( min1.first != v1_in_neighbors_end && min2.first != v2_in_neighbors_end ) {
		if( g.target_of( min1.first ) < g.target_of( min2.first ) ) {
			if( g.target_of( min1.first ) != v2 ) {
				++temp_edge_types.insert( std::make_pair( min1.second.first < min1.second.second ? std::make_pair( min1.second.first, min1.second.second ) : std::make_pair( min1.second.second, min1.second.first ), 0 ) ).first->second;
				++gaps;
				v3Vertices_end->first = g.target_of( min1.first );
				v3Vertices_end->second = connectivity;
				v3Vertices_end->second( 0, 2 ) = min1.second.first;
				v3Vertices_end->second( 2, 0 ) = min1.second.second;
				++v3Vertices_end;
			}
			min1 = next_union_element( v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it, v1_in_neighbors_end );
		} else if( g.target_of( min1.first ) > g.target_of( min2.first ) ) {
			if( g.target_of( min2.first ) != v1 ) {
				++temp_edge_types.insert( std::make_pair( min2.second.first < min2.second.second ? std::make_pair( min2.second.first, min2.second.second ) : std::make_pair( min2.second.second, min2.second.first ), 0 ) ).first->second;
				++gaps;
				v3Vertices_end->first = g.target_of( min2.first );
				v3Vertices_end->second = connectivity;
				v3Vertices_end->second( 1, 2 ) = min2.second.first;
				v3Vertices_end->second( 2, 1 ) = min2.second.second;
				++v3Vertices_end;
			}
			min2 = next_union_element( v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it, v2_in_neighbors_end );
		} else { // the next neighbor is shared by both v1 and v2, so it cannot be either and we do not need to check to exclude it
			++temp_edge_types.insert( std::make_pair( min1.second.first < min1.second.second ? std::make_pair( min1.second.first, min1.second.second ) : std::make_pair( min1.second.second, min1.second.first ), 0 ) ).first->second;
			++temp_edge_types.insert( std::make_pair( min2.second.first < min2.second.second ? std::make_pair( min2.second.first, min2.second.second ) : std::make_pair( min2.second.second, min2.second.first ), 0 ) ).first->second;
		 	v3Vertices_end->first = g.target_of( min1.first );
		 	v3Vertices_end->second = connectivity;
		 	v3Vertices_end->second( 0, 2 ) = min1.second.first;
		 	v3Vertices_end->second( 2, 0 ) = min1.second.second;
		 	v3Vertices_end->second( 1, 2 ) = min2.second.first;
		 	v3Vertices_end->second( 2, 1 ) = min2.second.second;
			++v3Vertices_end;
			min1 = next_union_element( v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it, v1_in_neighbors_end );
			min2 = next_union_element( v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it, v2_in_neighbors_end );
		}
	}
	while( min1.first != v1_in_neighbors_end ) {
		if( g.target_of( min1.first ) != v2 ) {
			++temp_edge_types.insert( std::make_pair( min1.second.first < min1.second.second ? std::make_pair( min1.second.first, min1.second.second ) : std::make_pair( min1.second.second, min1.second.first ), 0 ) ).first->second;
			++gaps;
			v3Vertices_end->first = g.target_of( min1.first );
			v3Vertices_end->second = connectivity;
			v3Vertices_end->second( 0, 2 ) = min1.second.first;
			v3Vertices_end->second( 2, 0 ) = min1.second.second;
			++v3Vertices_end;
		}
		min1 = next_union_element( v1_out_neighbors_it, v1_out_neighbors_end, v1_in_neighbors_it, v1_in_neighbors_end );
	}
	while( min2.first != v2_in_neighbors_end ) {
		if( g.target_of( min2.first ) != v1 ) {
			++temp_edge_types.insert( std::make_pair( min2.second.first < min2.second.second ? std::make_pair( min2.second.first, min2.second.second ) : std::make_pair( min2.second.second, min2.second.first ), 0 ) ).first->second;
			++gaps;
			v3Vertices_end->first = g.target_of( min2.first );
			v3Vertices_end->second = connectivity;
			v3Vertices_end->second( 1, 2 ) = min2.second.first;
			v3Vertices_end->second( 2, 1 ) = min2.second.second;
			++v3Vertices_end;
		}
		min2 = next_union_element( v2_out_neighbors_it, v2_out_neighbors_end, v2_in_neighbors_it, v2_in_neighbors_end );
	}

	std::size_t v3_count( v3Vertices_end-v3Vertices_begin );
	std::size_t v4_count( 0 );
	for( std::pair<const_vertex_iterator,connectivity_matrix >* it1( v3Vertices_begin ); it1 != v3Vertices_end; ++it1 ) { // for each v3 vertex computed above
		const_edge_iterator v3_out_neighbors_it( g.out_neighbors_begin( it1->first ) );
		const_edge_iterator v3_out_neighbors_end( g.out_neighbors_end( it1->first ) );
		const_edge_iterator v3_in_neighbors_it( g.in_neighbors_begin( it1->first ) );
		const_edge_iterator v3_in_neighbors_end( g.in_neighbors_end( it1->first ) );
		unsigned long v4_local_count( 0 ); // keep track of how many v4 vertices are only the result of the neighbors of this v3
		std::pair<const_edge_iterator,std::pair<connectivity_address_type,connectivity_address_type> > min( next_union_element( v3_out_neighbors_it, v3_out_neighbors_end, v3_in_neighbors_it, v3_in_neighbors_end ) );
		for( std::pair<const_vertex_iterator,connectivity_matrix>* it2( v3Vertices_begin ); it2 != v3Vertices_end; ++it2 ) {	
			while( min.first != v3_in_neighbors_end && g.target_of( min.first ) < it2->first ) {
				if( g.target_of( min.first ) != v1 && g.target_of( min.first ) != v2 ) {
					++temp_edge_types.insert( std::make_pair( min.second.first < min.second.second ? std::make_pair( min.second.first, min.second.second ) : std::make_pair( min.second.second, min.second.first ), 0 ) ).first->second;
					++v4_local_count;
					it1->second( 0, 3 ) = 0;
					it1->second( 3, 0 ) = 0;
					it1->second( 1, 3 ) = 0;
					it1->second( 3, 1 ) = 0;
					it1->second( 2, 3 ) = min.second.first;
					it1->second( 3, 2 ) = min.second.second;
					++counts.insert( std::make_pair( mapper.canonical_subgraph_address( it1->second ), 0 ) ).first->second;
				}
				min = next_union_element( v3_out_neighbors_it, v3_out_neighbors_end, v3_in_neighbors_it, v3_in_neighbors_end );
			}
			if( min.first == v3_in_neighbors_end || g.target_of( min.first ) > it2->first ) {
				if( it1->first < it2->first ) {
					++gaps;
					it1->second( 0, 3 ) = it2->second( 0, 2 );
					it1->second( 3, 0 ) = it2->second( 2, 0 );
					it1->second( 1, 3 ) = it2->second( 1, 2 );
					it1->second( 3, 1 ) = it2->second( 2, 1 );
					it1->second( 2, 3 ) = 0;
					it1->second( 3, 2 ) = 0;
					++counts.insert( std::make_pair( mapper.canonical_subgraph_address( it1->second ), 0 ) ).first->second;
				}
			} else {
				if( it1->first < it2->first ) {
					++temp_edge_types.insert( std::make_pair( min.second.first < min.second.second ? std::make_pair( min.second.first, min.second.second ) : std::make_pair( min.second.second, min.second.first ), 0 ) ).first->second;
					it1->second( 0, 3 ) = it2->second( 0, 2 );
					it1->second( 3, 0 ) = it2->second( 2, 0 );
					it1->second( 1, 3 ) = it2->second( 1, 2 );
					it1->second( 3, 1 ) = it2->second( 2, 1 );
					it1->second( 2, 3 ) = min.second.first;
					it1->second( 3, 3 ) = min.second.second;
					++counts.insert( std::make_pair( mapper.canonical_subgraph_address( it1->second ), 0 ) ).first->second;
				}
				min = next_union_element( v3_out_neighbors_it, v3_out_neighbors_end, v3_in_neighbors_it, v3_in_neighbors_end );
			}
		}
		while( min.first != v3_in_neighbors_end ) {
			if( g.target_of( min.first ) != v1 && g.target_of( min.first ) != v2 ) {
				++temp_edge_types.insert( std::make_pair( min.second.first < min.second.second ? std::make_pair( min.second.first, min.second.second ) : std::make_pair( min.second.second, min.second.first ), 0 ) ).first->second;
				++v4_local_count;
				it1->second( 0, 3 ) = 0;
				it1->second( 3, 0 ) = 0;
				it1->second( 1, 3 ) = 0;
				it1->second( 3, 1 ) = 0;
				it1->second( 2, 3 ) = min.second.first;
				it1->second( 3, 2 ) = min.second.second;
				++counts.insert( std::make_pair( mapper.canonical_subgraph_address( it1->second ), 0 ) ).first->second;
			}
			min = next_union_element( v3_out_neighbors_it, v3_out_neighbors_end, v3_in_neighbors_it, v3_in_neighbors_end );
		}
		v4_count += v4_local_count;
		gaps += 2*v4_local_count;
		it1->second( 0, 3 ) = 0;
		it1->second( 3, 0 ) = 0;
		it1->second( 1, 3 ) = 0;
		it1->second( 3, 1 ) = 0;
		it1->second( 2, 3 ) = 0;
		it1->second( 3, 2 ) = 0;
		counts.insert( std::make_pair( mapper.canonical_subgraph_address( it1->second ), 0 ) ).first->second += g.vertex_count() - 2 - ((v3Vertices_end-v3Vertices_begin) + v4_local_count);
	}
		
	for( typename std::map<std::pair<connectivity_address_type,connectivity_address_type>,unsigned long>::const_iterator it( edge_types.begin() ); it != edge_types.end(); ++it ) {
		connectivity( 2, 3 ) = it->first.first; // THESE ARE ALWAYS GOING TO BE COLLAPSED ISOMORPHICALLY EQUIVALENTLY HERE, SO WE CAN BENCHMARK AFTER GETTING RID OF THE UGLY MAKE_PAIRS
		connectivity( 3, 2 ) = it->first.second;
		typename std::map<std::pair<connectivity_address_type,connectivity_address_type>,unsigned long>::const_iterator temp_it( temp_edge_types.find( it->first ) );
		unsigned long count( it->second );
		if( temp_it != temp_edge_types.end() ) {
			count -= temp_it->second;
			if( it->first.first == 0 && it->first.second == 0 ) {
				count -= !static_cast<bool>( connectivity( 0, 1 ) + connectivity( 1, 0 ) ) + (2 + v3_count) * (g.vertex_count() - 2 - v3_count) - 3 * v4_count;
			}
		}
		counts.insert( std::make_pair( mapper.canonical_subgraph_address( connectivity ), 0 ) ).first->second += count;
	}

	
	return counts;
}

}

#endif
