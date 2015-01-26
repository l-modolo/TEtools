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

import argparse
import configparser
from os import path
from os import stat
from os import mkdir
from subprocess import Popen
from subprocess import PIPE
from select import select
import atexit

tools_path = path.realpath(__file__)[:-len('PingPong.py')]
config = configparser.ConfigParser()
if not path.isfile(tools_path+'PingPong.ini'):
    print("'PingPong.ini' file not found, writing default one.")
    config['programs'] = {'pingpongpro' : tools_path+'pingpongpro/pingpongpr',
                          'bowtie' : tools_path+'bowtie/bowtie',
                          'min-alignment-length' : '10',
                          'max-alignment-length' : '40',
                          'TE_file' : 'pingpong_working_directory/TE.fasta',
                          'csv_file' : 'pingpong_working_directory/TE.csv',
                          'sam_file' : 'pingpong_working_directory/TE.sam',
                          'output_folder' : 'pingpong_working_directory/',
                          'thread' : '4'
                         }
    with open(tools_path+'PingPong.ini', 'w') as configfile:
        config.write(configfile)
config.read(tools_path+'PingPong.ini')

parser = argparse.ArgumentParser(prog='PingPong.py')
parser.add_argument('-fastqs', action='store', dest='fastq_files', help='RNA sam file(s)', nargs='*')
parser.add_argument('-fasta', action='store', dest='fasta_file', default = None, help='TE sequence fasta file')
parser.add_argument('-TE_names', action='store', dest='TE_names', default=500, help='names of the TE family to check', nargs='*')
parser.add_argument('-version', action='store_true', default=False, dest='version', help='print version')
args = parser.parse_args()

if args.version:
    print("1.0.0")
    exit(0)

class Pingpong:
    procs = list()

    # we want to kill every subprocess if we kill countTE.py
    def kill_subprocesses(self):
        if len(self.procs) > 0:
            for proc in self.procs:
                proc.kill()

    @staticmethod
    def __output_reader(proc, retVal=''): 
        while (select([proc.stdout],[],[],0)[0]!=[]):   
            retVal+=str(proc.stdout.read(1))
        return retVal

    def __init__(self, config, fasta_file, TE_names, fastq_files):
        self.fasta_file = str(fasta_file)
        if not path.isfile(self.fasta_file):
            print("error: cannot open fasta file")
            exit(1)
        self.TE_by_family = True
        if len(TE_names) >= 1:
            self.TE_nams = TE_names
        else:
            print("error: need TE name")
        self.config = config
        self.fastq_files = fastq_files
        self.TE_names = list()
        self.TE_sequences = list()
        self.__check_TE()
        self.__write_csv()
        self.__mapping_index()
        self.__mapping()

    def __check_TE(self):
        with open(self.fasta_file) as fasta_file_handle:
            is_found = False
            TE_sequence = str()
            for line in fasta_file_handle:
                if is_found:
                    self.TE_sequences[len(self.TE_sequences)-1] += line
                if line[0] == '>':
                    if any(line.find(TE_name) != -1 for TE_name in self.TE_names):
                        is_found = True
                        self.TE_names.append(line[1:])
                        self.TE_sequences.append('')
                    else:
                        is_found = False
            print("number of copies found: "+str(len(self.TE_names)))
            with open(self.config['TE_file'], "w") as TE_file_handle:
                TE_file_handle.write(TE_sequence)

    def __write_csv(self):
        with open(self.config['csv_file'], "w") as csv_file_handle:
            for i in range(len(self.TE_names)):
                csv_file_handle.write(str(self.TE_names[i])+ ',+,'+ str(self.TE_names[i])+ ',0,'+ str(len(self.TE_sequences[i]))+ "\n")

    def __mapping_index(self):
        print('building index')
        if path.isfile(self.config['TE_file']):
            self.index_file = self.config['TE_file']+'.index'
            bowtie_cmd = list()
            bowtie_cmd.append(self.config['bowtie']+'-build')
            bowtie_cmd += ['-f', self.config['TE_file'], self.index_file]
            print(bowtie_cmd)
            self.procs.append(Popen(bowtie_cmd, stdout=PIPE, stderr=PIPE))
            self.procs[len(self.procs)-1].wait()
            if self.procs[len(self.procs)-1].returncode != 0:
                for line in self.procs[len(self.procs)-1].stderr:
                    print(str(line))
                self.procs.pop()
                exit(1)
            for line in self.procs[len(self.procs)-1].stdout:
                print(line)
            self.procs.pop()
        else:
            print(str(self.config['TE_file'])+' file not found')
            exit(1)

    def __mapping(self):
        print('mapping reads')
        bowtie_cmd = list()
        bowtie_cmd.append(str(self.config['bowtie']))
        bowtie_cmd += ['-S',
                        '-p '+str(self.config['thread']),
                        '--time',
                        '--best',
                        str(self.index_file)]
        fastq_list = '-s '
        for fastq_file in self.fastq_files:
            if path.isfile(fastq_file):
                str(self.fastq_file)+','
            else:
                print('error: fastq file '+str(fastq_file)+' not found')
                exit(1)
        bowtie_cmd += [fastq_list[:-1], str(self.config['sam_file'])]
        print(bowtie_cmd)
        self.procs.append(Popen(bowtie_cmd, stdout=PIPE, stderr=PIPE))
        self.procs[len(self.procs)-1].wait()
        if self.procs[len(self.procs)-1].returncode != 0:
            for line in self.procs[len(self.procs)-1].stderr:
                print(line)
            self.procs.pop()
            exit(1)
        for line in self.procs[len(self.procs)-1].stdout:
            print(line)
        self.procs.pop()

    def __pingpong(self):
        print('computing pingpong metrics')
        pingpongpro_cmd = list()
        pingpongpro_cmd.append(str(self.config['pingpongpro']))
        pingpongpro_cmd += ['-b',
                        '-i '+str(self.config['sam_file']),
                        '-l '+str(self.config['min-alignment-length']),
                        '-L '+str(self.config['max-alignment-length']),
                        '-m weighted ',
                        '-o '+str(self.config['output_folder']),
                        '-t '+str(self.config['csv_file']),
                        '-v']
        pingpongpro_cmd += [fastq_list[:-1], str(self.config['sam_file'])]
        print(pingpongpro_cmd)
        self.procs.append(Popen(pingpongpro_cmd, stdout=PIPE, stderr=PIPE))
        self.procs[len(self.procs)-1].wait()
        if self.procs[len(self.procs)-1].returncode != 0:
            for line in self.procs[len(self.procs)-1].stderr:
                print(line)
            self.procs.pop()
            exit(1)
        for line in self.procs[len(self.procs)-1].stdout:
            print(line)
        self.procs.pop()

@atexit.register
def kill_subprocesses():
    Pingpong.kill_subprocesses()

TE = Pingpong(config['programs'], args.fasta_file, args.TE_names, args.fastq_files)
