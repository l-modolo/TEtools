#!/usr/bin/python3
# -*-coding:Utf-8 -*

#Copyright (C) 2013 Laurent Modolo

#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU Lesser General Public
#License as published by the Free Software Foundation; either
#version 2.1 of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#Lesser General Public License for more details.

#You should have received a copy of the GNU Lesser General Public
#License along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

import configparser
from os import path
import argparse
import subprocess

tools_path = path.realpath(__file__)[:-7]

config = configparser.ConfigParser()
if not path.isfile(tools_path+'UrQt.ini'):
    print("'UrQt.ini' file not found, writing default one.")
    config['programs'] = {'urqt' : tools_path+'UrQt',
                          'thread' : '4'
                         }
    with open(tools_path+'UrQt.ini', 'w') as configfile:
        config.write(configfile)
config.read(tools_path+'UrQt.ini')

parser = argparse.ArgumentParser(prog='UrQt.py')
parser.add_argument('-in', action='store', dest='input_file', help='fastq input file')
parser.add_argument('-out', action='store', dest='output_file', help='fastq output file')
parser.add_argument('-inpair', action='store', dest='input_file_pair', default=None, help='fastq input pair file')
parser.add_argument('-outpair', action='store', dest='output_file_pair',default=None,  help='fastq output pair file')
parser.add_argument('-t', action='store', default="20", dest='threshold', help='threshold phred')
parser.add_argument('-pos', action='store', default="both", dest='pos', help='expected position of trimmed sequence in the read (default: both)')
parser.add_argument('-min_read_size', action='store', dest='min_read_size', default=None, help='remove all reads smaller than this size after the trimming step (default: 0)')
parser.add_argument('-polyN', action='store', dest='polyN', default=None, help='polyN to trim (default: QC trimming)')
parser.add_argument('-r', action='store', dest='empty_read', default=None, help='polyN to trim (default: QC trimming)')
parser.add_argument('-version', action='store_true', dest='version', default=False, help='version')
args = parser.parse_args()

if args.version:
    print("1.0.18")
    exit(0)

input_file = str(args.input_file)#[2:]
# input_file = input_file[:(len(input_file)-2)]

urqt_cmd = str(config['programs']['urqt'])+' --v --m '+str(config['programs']['thread'])+' --t '+str(args.threshold)+' --pos '+str(args.pos)+' --in '+str(input_file)+' --out '+str(args.output_file)

if args.input_file_pair is not None and args.output_file_pair is not None:
    input_file_pair = str(args.input_file_pair)[2:]
    input_file_pair = input_file_pair[:(len(input_file_pair)-2)]
    urqt_cmd += ' --inpair '+str(input_file_pair)+' --outpair '+str(args.output_file_pair)

if args.min_read_size is not None:
    urqt_cmd += ' --min_read_size '+str(args.min_read_size)

if args.polyN is not None:
    urqt_cmd += ' --polyN '+str(args.polyN)

print(urqt_cmd)
urqt_process = subprocess.Popen(str(urqt_cmd), shell=True)
urqt_process.wait()
exit(0)