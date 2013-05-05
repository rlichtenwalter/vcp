/*
Copyright (C) 2013 by Ryan N. Lichtenwalter
Email: rlichtenwalter@gmail.com

This file is part of the Vertex Collocation Profiles code base.

The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VCP_SQUARE_MATRIX_H
#define VCP_SQUARE_MATRIX_H

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace vcp {

template <typename value_type,std::size_t n = 0>
class square_matrix {
	public:
		std::size_t size() const;
		value_type & operator()( std::size_t row, std::size_t column );
		value_type const & operator()( std::size_t row, std::size_t column ) const;
		template <typename value_type_,std::size_t n_> friend std::ostream & operator<<( std::ostream & os, square_matrix<value_type_,n_> const & matrix );
	private:
		std::array<value_type,n*n> data = {{0}};
};

template <typename value_type,size_t n> std::ostream & operator<<( std::ostream & os, square_matrix<value_type,n> const & matrix );

template <typename value_type,std::size_t n>
std::size_t square_matrix<value_type,n>::size() const {
	return std::sqrt( data.size() );
}

template <typename value_type,std::size_t n>
value_type & square_matrix<value_type,n>::operator()( std::size_t row, std::size_t column ) {
	return const_cast<value_type &>( static_cast<square_matrix<value_type,n> const &>(*this)( row, column ) );
}

template <typename value_type,std::size_t n>
value_type const & square_matrix<value_type,n>::operator()( std::size_t row, std::size_t column ) const {
	return data[size()*row + column];
}

template <typename value_type,std::size_t n>
std::ostream & operator<<( std::ostream & os, square_matrix<value_type,n> const & matrix ) {
	for( std::size_t row( 0 ); row < matrix.size(); ++row ) {
		for( std::size_t column( 0 ); column < matrix.size(); ++column ) {
			os << matrix( row, column );
			if( column < matrix.size() - 1 ) {
				os << ',';
			}
		}
		os << '\n';
	}
	return os;
}

template <typename value_type> std::ostream & operator<<( std::ostream & os, square_matrix<value_type,0> const & matrix );

template <typename value_type>
class square_matrix<value_type,0> {
	public:
		std::size_t size() const;
		void resize( std::size_t n );
		square_matrix() = default;
		square_matrix( std::size_t n );
		template <std::size_t n> square_matrix( square_matrix<value_type,n> const & matrix );
		template <std::size_t n> square_matrix<value_type,0> & operator=( square_matrix<value_type,n> const & matrix );
		value_type & operator()( std::size_t row, std::size_t column );
		value_type const & operator()( std::size_t row, std::size_t column ) const;
	private:
		std::vector<value_type> data;
};

template <typename value_type>
std::size_t square_matrix<value_type,0>::size() const {
	return std::sqrt( data.size() );
}

template <typename value_type>
void square_matrix<value_type,0>::resize( std::size_t n ) {
	data.resize( n * n );
}

template <typename value_type>
square_matrix<value_type,0>::square_matrix( std::size_t n ) : data( n*n, 0 ) {
}

template <typename value_type>
template <typename std::size_t n>
square_matrix<value_type,0>::square_matrix( square_matrix<value_type,n> const & matrix ) : data( n*n ) {
	std::copy( &matrix( 0, 0 ), &matrix( n - 1, n - 1 ), &data[0] );
}

template <typename value_type>
template <typename std::size_t n>
square_matrix<value_type,0> & square_matrix<value_type,0>::operator=( square_matrix<value_type,n> const & matrix ) {
	if( this != &matrix ) {
		data.resize( n );
		std::copy( &matrix.data[0], &matrix.data[n*n], &data[0] );
	}
	return *this;
}

template <typename value_type>
value_type & square_matrix<value_type,0>::operator()( std::size_t row, std::size_t column ) {
	return const_cast<value_type &>( static_cast<square_matrix<value_type,0> const &>(*this)( row, column ) );
}

template <typename value_type>
value_type const & square_matrix<value_type,0>::operator()( std::size_t row, std::size_t column ) const {
	return data[size()*row + column];
}

template <typename value_type>
std::ostream & operator<<( std::ostream & os, square_matrix<value_type,0> const & matrix ) {
	for( std::size_t row( 0 ); row < matrix.size(); ++row ) {
		for( std::size_t column( 0 ); column < matrix.size(); ++column ) {
			os << matrix( row, column );
			if( column < matrix.size() - 1 ) {
				os << ',';
			}
		}
		os << '\n';
	}
	return os;
}

}

#endif
