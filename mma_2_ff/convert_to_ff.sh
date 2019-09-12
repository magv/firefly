#!/bin/bash
#==================================================================================
#    FireFly - Reconstructing rational functions and polynomial over finite fields.
#    Copyright (C) 2019  Jonas Klappert and Fabian Lange
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#==================================================================================

# @Obsolete
# Script to convert a mathematica list to compilabale C++ code
mkdir -p ff_conv

funs=$1
vars=$2
threads=$3

math -run "functions=\"$funs\";variables=\"$vars\";nthreads=$threads" -script convert_to_ff.m

for f in ff_conv/*.cpp; do
    echo "Converting $f"
    sed -i 's/"mpz_class(",/mpz_class(/g' $f
    sed -i 's/,")"/)/g' $f
done

