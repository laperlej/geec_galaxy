import sys
import argparse
import subprocess
import tempfile
import itertools
import os
import json
import os.path
import h5py
import scipy.stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'geec_tools/python/utils'))
import config

PUBLIC_DATA_ROOT = "/home/laperlej/geec/public"

class Wig(object):
    def __init__(self, wigfile):
        self.file = open(wigfile, 'r')
        self.chromsize = {}
        self.cursor = {}
        self.reset_cursor()

    def __str__(self):
        string = []
        for chrom, size in self.chromsize.iteritems():
            string.append("{0} {1}".format(chrom, size))
        return "\n".join(string)

    def reset_cursor(self):
        self.cursor["type"] = None
        self.cursor["count"] = 0
        self.cursor["position"] = 0
        self.cursor["chrom"] = None
        self.cursor["start"] = 0
        self.cursor["step"] = 1
        self.cursor["span"] = 1

    def read(self):
        for line in self.file:
            line = line.strip()
            if line[0].isdigit():
                self.cursor["position"] = int(line.split()[0])
                self.cursor["count"] += 1
            elif line.startswith("fixedStep"):
                if self.cursor["type"] is not None:
                    self.update_chromsize()
                self.cursor["type"] = "fixed"
                self.parse_attributes(line)
            elif line.startswith("variableStep"):
                if self.cursor["type"] is not None:
                    self.update_chromsize()
                self.cursor["type"] = "variable"
                self.parse_attributes(line)
        if self.cursor["type"] is not None:
            self.update_chromsize()


    def parse_attribute(self, tag):
        return tag.split("=")[1].strip()

    def parse_attributes(self, line):
        line = line.strip().split()
        for tag in line:
            if tag.startswith("chrom"):
                self.cursor["chrom"] = self.parse_attribute(tag)
            elif tag.startswith("start"):
                self.cursor["start"] = int(self.parse_attribute(tag))
            elif tag.startswith("step"):
                self.cursor["step"] = int(self.parse_attribute(tag))
            elif tag.startswith("span"):
                self.cursor["span"] = int(self.parse_attribute(tag))

    def update_chromsize(self):
        size = 0
        if self.cursor["type"] == "fixed":
            size = self.cursor["start"] + self.cursor["count"] * self.cursor["step"] + self.cursor["span"]
        elif self.cursor["type"] == "variable":
            size = self.cursor["position"] + self.cursor["span"]
        if size > self.chromsize.get(self.cursor["chrom"],0):
            self.chromsize[self.cursor["chrom"]] = size
        self.reset_cursor()

def bw_to_hdf5(raw_file, name, assembly, user_hdf5, resolution):
    """Usage: to_hdf5 {dataset.bw}
                      {name}
                      {chrom_sizes}
                      {output.hdf5}
                      {bin_size}\n"""
    arguments = [config.BW_TO_HDF5,
                 raw_file,
                 name,
                 config.get_chrom_sizes(assembly),
                 user_hdf5,
                 resolution]
    subprocess.call(arguments)

def bg_to_hdf5(raw_file, name, assembly, user_hdf5, resolution):
    """Usage: to_hdf5 {dataset.bw}
                      {name}
                      {chrom_sizes}
                      {output.hdf5}
                      {bin_size}\n"""
    arguments = [config.BG_TO_HDF5,
                 raw_file,
                 name,
                 config.get_chrom_sizes(assembly),
                 user_hdf5,
                 resolution]
    subprocess.call(arguments)

def wig_to_bigwig(wig_file, bigwig_file):
    # wigToBigWig in.wig chrom.sizes out.bw
    wig = Wig(wig_file)
    wig.read()
    chromsizes_file = temp_name()
    with open(chromsizes_file, "w") as chrom_size:
        chrom_size.write(str(wig))
    arguments = [config.WIG_TO_BW,
                 wig_file,
                 chromsizes_file,
                 bigwig_file]
    subprocess.call(arguments)


def filter_hdf5(name, assembly, user_hdf5, filtered_hdf5, resolution, include, exclude):
    """Usage: filter    {input.hdf5}
                        {name}
                        {output.hdf5}
                        {chrom_sizes}
                        {bin_size}
                        {include.bed}
                        {exclude.bed}\n");"""
    arguments = [config.FILTER,
                 user_hdf5,
                 name,
                 filtered_hdf5,
                 config.get_chrom_sizes(assembly),
                 resolution,
                 include,
                 exclude]
    subprocess.call(arguments)

def correlate(input_list, assembly, correlation_file, resolution):
    """Usage: correlation {input_list}
                          {chrom_sizes}
                          {output.results}
                          {bin_size}\n");"""
    arguments = [config.CORRELATION,
                 input_list,
                 config.get_chrom_sizes(assembly),
                 correlation_file,
                 resolution]
    subprocess.call(arguments)

def can_slice_matrix(md5s, files, assembly, resolution, include, exclude, metric):
    # Slice if md5s and not files and is Pearson
    # if .mat exist
    if md5s and not files and metric == "pearson":
        path = config.get_matrix(assembly, resolution, include, exclude) 
        if os.path.isfile(path):
            return True
    return False

def is_nm(md5s, files):
    # verify if nm
    return bool(md5s and files and metric == "pearson")

def correlate_nm(input_list1, input_list2, assembly, correlation_file, resolution):
    """Usage: correlation_nm {input_list1}
                             {input_list2}
                             {chrom_sizes}
                             {output.results}
                             {bin_size}\n");"""
    subprocess.call([config.CORRELATION_NM,
                     input_list1,
                     input_list2,
                     config.get_chrom_sizes(assembly),
                     correlation_file,
                     resolution
                     ])

def make_matrix_nm(input_list1, input_list2, correlation_file, precalc_matrix, output_matrix, meta_json = ""):
    """
    python make_matrix.py {list_path} {chrom_size} {corr_path} {output_path}
    """
    arguments = ['python',
                 config.MAKE_MATRIX_NM,
                 input_list1,
                 input_list2,
                 correlation_file,
                 precalc_matrix,
                 output_matrix
                 ]
    if meta_json:
        arguments += [meta_json]
    subprocess.call(arguments)

def slice_matrix(md5s, assembly, resolution, include, exclude, output):
    """
    python geec_slice_md5sum.py matrix.mat md5_1 md5_2 md5_N > output.mat
    """
    arguments = ['python',
                 config.GEEC_SLICE_MD5SUM,
                 config.get_matrix(assembly, resolution, include, exclude)]
    arguments += md5s

    with open(output, 'w') as output_file:
        subprocess.call(arguments, stdout=output_file)

def make_matrix(input_list, correlation_file, output_matrix, meta_json = ""):
    """
    python make_matrix.py {list_path} {chrom_size} {corr_path} {output_path}
    """
    arguments = ['python', 
                 config.MAKE_MATRIX,
                 input_list,
                 correlation_file,
                 output_matrix]
    if meta_json:
      arguments += [meta_json]
    subprocess.call(arguments)

def rank_hdf5(hdf5_path):
  h5f = h5py.File(hdf5_path, 'r+')
  for group in h5f:
    md5 = group
  for dset in h5f[md5]:
    data = h5f[md5][dset]
    assert(n < 3024616) #overflow on int64
    n = len(data)
    data[...] = scipy.stats.rankdata(data, method="ordinal")
    data.attrs['sumX'] = n * (n + 1) / 2 #sum of natural numbers
    data.attrs['sumXX'] = n * (n + 1) * (2 * n + 1) / 6 #sum of perfect squares
  h5f.close()

def create_input_list(input_list):
    file_path = tmp_name()
    with open(file_path, 'w') as input_file:
        for path, label in input_list:
            input_file.write('{0}\t{1}\n'.format(path.strip(),label.strip()))
    return file_path

def parse_md5s(md5_json):
    md5s = []
    if md5_json:
      for md5 in md5_json["datasets"]:
        md5s.append(md5["md5sum"])
    return md5s

def tmp_name():
    fd, temp_path = tempfile.mkstemp()
    os.close(fd)
    os.remove(temp_path)
    return temp_path

def listjson2dictjson(old_json):
    new_json = {"datasets":{}}
    for token in old_json.get("datasets", []):
        new_json["datasets"][token["md5sum"]] = token
    return new_json

def main():
    """
    sample input:
    Namespace(bigwigs=['/Users/Jon/Work/galaxy-dist/database/files/000/dataset_172.dat',
                       '/Users/Jon/Work/galaxy-dist/database/files/000/dataset_173.dat'],
              bin=10000,
              exclude='merged',
              include='none',
              labels=['MS000501.cd4phtc.H3K4me3.signal.bigWig',
                      'MS010002.monocyte.H3K27me3.signal.bigWig'],
              md5s='/Users/Jon/Work/galaxy-dist/database/files/000/dataset_169.dat',
              output='/Users/Jon/Work/galaxy-dist/database/files/000/dataset_246.dat')
    """
    #parse arguments
    parser = argparse.ArgumentParser(description='GeEC interface for galaxy')
    parser.add_argument('--files', nargs='*')
    parser.add_argument('--types', nargs='*')
    parser.add_argument('--labels', nargs='*')
    parser.add_argument('--md5s')
    parser.add_argument('--include')
    parser.add_argument('--exclude')
    parser.add_argument('--bin')
    parser.add_argument('--output')
    parser.add_argument('--assembly')
    parser.add_argument('--metric')
    args = parser.parse_args()
    if args.md5s == "None":
        args.md5s = []
    if args.files == ['None']:
        args.files = []
    if args.types == ['None']:
        args.files = []
    if args.labels == ['None']:
        args.labels = []

    #read md5 list file for public data
    if args.md5s:
      md5_json = json.load(open(args.md5s))
      md5s = parse_md5s(md5_json)
    else:
      md5_json = {}
      md5s = []

    md5_json = listjson2dictjson(md5_json)

    #slice_matrix or make_matrix
    if can_slice_matrix(md5s, args.files, args.assembly, args.bin, args.include, args.exclude, args.metric):
        slice_matrix(md5s, args.assembly, args.bin, args.include, args.exclude, args.output)
    else:
        #public data paths
        include_path = config.get_region(args.assembly, args.include)
        exclude_path = config.get_region(args.assembly, args.exclude)
        #create temporary input files for geec executables
        input_list1 = []
        input_list2 = []

        user_input_list = []
        for file, datatype, label in itertools.izip(args.files, args.types, args.labels):
            user_hdf5 = tmp_name()
            user_filtered_hdf5 = tmp_name()
            label = label.split("/")[-1]
            user_input_list.append((file, datatype, label, user_hdf5, user_filtered_hdf5))
            input_list1.append((user_filtered_hdf5, label))

        for md5 in md5s:
            hdf5_path = config.get_hdf5(md5, args.assembly, args.bin, args.include, args.exclude, args.metric)
            if os.path.isfile(hdf5_path):
                input_list2.append((hdf5_path, md5))
            else:
                print "{0} is missing".format(md5_json["datasets"][md5].get("file_name", "unknown"))

        correlation_file = tmp_name()

        # convert user bigwigs to hdf5 and filter it
        for raw_file, datatype, name, user_hdf5, user_filtered_hdf5 in user_input_list:
            if datatype.lower() == "bigwig":
                bw_to_hdf5(raw_file, name, args.assembly, user_hdf5, args.bin)
            elif datatype.lower() == "bedgraph":
                bg_to_hdf5(raw_file, name, args.assembly, user_hdf5, args.bin)
            elif datatype.lower() == "wig":
                tmp_file = tmp_name()
                wig_to_bigwig(raw_file, tmp_file)
                bw_to_hdf5(tmp_file, name, args.assembly, user_hdf5, args.bin)
            else:
                print "Could not determine type for {0}".format(name)
                continue
            filter_hdf5(name, args.assembly, user_hdf5, user_filtered_hdf5, args.bin, include_path, exclude_path)
            if args.metric == "spearman":
              rank_hdf5(user_filtered_hdf5)

        if is_nm(md5s, args.files):
            input_list_path1 = create_input_list(input_list1)
            input_list_path2 = create_input_list(input_list2)
            #correlate all uncorrelated matrix cells
            correlate_nm(input_list_path1, input_list_path2, args.assembly, correlation_file, args.bin)

            #generate the final matrix
            precalc_matrix = config.get_matrix(args.assembly, args.bin, include_path, exclude_path)
            make_matrix_nm(input_list_path1, input_list_path2, correlation_file, precalc_matrix, args.output, args.md5s)
        else:
            input_list_path = create_input_list(input_list1 + input_list2)
            #correlate all uncorrelated matrix cells
            correlate(input_list_path, args.assembly, correlation_file, args.bin)

            #generate the final matrix
            make_matrix(input_list_path, correlation_file, args.output, args.md5s)

if __name__ == '__main__':
    main()
