#!/usr/bin/python3
# -*-coding:Utf-8 -*

# Copyright (C) 2015 Laurent Modolo
 
# This file is part of TEtools suite for galaxy.

# TEtools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# TEtools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with TEtools.  If not, see <http://www.gnu.org/licenses/>.

import argparse

parser = argparse.ArgumentParser(prog='countTE.py')
parser.add_argument('-rosette', action='store', dest='rosette_file', help='Rosette file to fix')
parser.add_argument('-column', action='store', dest='fix_column', default=2, help='Rosette column to fix')
parser.add_argument('-fix', action='store', dest='fix_file', help='Rosette column to group count')
parser.add_argument('-output', action='store', dest='output_file', help='Rosette file fixed')

args = parser.parse_args()

fix_dict = dict()
with open(str(args.fix_file)) as fix_file_handle:
    for line in fix_file_handle:
        line = line.split()
        if len(line) == 2:
            fix_dict[line[0]] = line[1]

fix_column = int(args.fix_column)-1
with open(str(args.rosette_file)) as rosette_file_handle, open(str(args.output_file), "w") as output_file_handle:
    for line in rosette_file_handle:
        line = line.split()
        if line[fix_column] in fix_dict:
            line[fix_column] = fix_dict[line[fix_column]]
        output_file_handle.write(' '.join(line)+"\n")
