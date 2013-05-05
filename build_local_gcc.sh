#! /bin/sh

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

if test $# -ne 2 ; then
	echo "usage: $0 <gcc_version (e.g. 4.8.0)> <target (e.g. $HOME/opt/)>" 1>&2
	exit 1
fi

version=$1
target=$2

cd $target || exit 1
wget http://www.netgull.com/gcc/releases/gcc-${version}/gcc-${version}.tar.gz || exit 1
tar -xzf gcc-${version}.tar.gz || exit 1
cd gcc-${version} || exit 1
./contrib/download_prerequisites || exit 1
cd .. || exit 1
mkdir objdir || exit 1
cd objdir || exit 1
$PWD/../gcc-${version}/configure --prefix=${target}/gcc-${version} || exit 1
make || exit 1
make install || exit 1

echo "You should add the following to your environment (e.g. bash shell, 64-bit system):"
echo "    export PATH=~/opt/gcc-${version}/bin:\$PATH"
echo "    export LD_LIBRARY_PATH=${target}/gcc-${version}/lib64:\$LD_LIBRARY_PATH"

