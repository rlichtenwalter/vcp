#! /usr/bin/perl

##
# Copyright (C) 2013 by Ryan N. Lichtenwalter
# Email: rlichtenwalter@gmail.com
#
# This file is part of the Vertex Collocation Profiles code base.
#
# The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
##

use strict;
use warnings FATAL => 'all';

my $edges_found = 0;
my $out_linenum = 0;
my $first_edge_of_line = 1;

my $line = <STDIN>;
$line =~ m/\*Vertices ([0-9]+)$/;
my $vertex_count = $1;

while( my $line = <STDIN> ) {
	if( !$edges_found && $line =~ m/^\*Edges.*$/ ) {
		$edges_found = 1;
		next;
	}
	if( $edges_found ) {
		$line =~ m/^([0-9]+) ([0-9]+) ([0-9]+)$/;
		while( $out_linenum < $1 ) {
			print "\n";
			++$out_linenum;
			$first_edge_of_line = 1;
		}
		if( !$first_edge_of_line ) {
			print " ";
		}
		print "$2,$3";
		$first_edge_of_line = 0;
	}
}

while( $out_linenum < $vertex_count ) {
	print "\n";
	++$out_linenum;
}

