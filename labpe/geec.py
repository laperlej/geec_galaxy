import sys
import argparse
import subprocess
import tempfile
import itertools
import os
import json

sys.path.insert(0, '/home/galaxy/geec_tools/python/utils')
import config

def to_hdf5(raw_file, name, assembly, user_hdf5, resolution):
    """Usage: to_hdf5 {dataset.bw}
                      {name}
                      {chrom_sizes}
                      {output.hdf5}
                      {bin_size}\n"""
    subprocess.call([config.TO_HDF5,
                     raw_file,
                     name,
                     config.CHROM_SIZE[assembly],
                     user_hdf5,
                     resolution])

def filter_hdf5(name, assembly, user_hdf5, filtered_hdf5, resolution, include, exclude):
    """Usage: filter    {input.hdf5}
                        {name}
                        {output.hdf5}
                        {chrom_sizes}
                        {bin_size}
                        {include.bed}
                        {exclude.bed}\n");"""
    subprocess.call([config.FILTER,
                     user_hdf5,
                     name,
                     filtered_hdf5,
                     config.CHROM_SIZE[assembly],
                     resolution,
                     include,
                     exclude
                    ])

def correlate(input_list, assembly, correlation_file, resolution):
    """Usage: correlation {input_list}
                          {chrom_sizes}
                          {output.results}
                          {bin_size}\n");"""
    subprocess.call([config.CORRELATION,
                     input_list,
                     config.CHROM_SIZE[assembly],
                     correlation_file,
                     resolution
                     ])

def make_matrix(input_list, assembly, correlation_file, output_matrix, meta_json = ""):
    """
    python make_matrix.py {list_path} {chrom_size} {corr_path} {output_path}
    """
    arguments = ['python', 
                 config.MAKE_MATRIX,
                 input_list,
                 config.CHROM_SIZE[assembly],
                 correlation_file,
                 output_matrix]
    if meta_json:
      arguments += meta_json
    subprocess.call(arguments)

def create_input_list(input_list):
    file_path = tmp_name()
    with open(file_path, 'w') as input_file:
        for path, label in input_list:
            input_file.write('{0}\t{1}\n'.format(path.strip(),label.strip()))
    return file_path

def parse_md5s(md5_json):
    md5s = []
    if md5_json:
      for md5 in md5_json["datasets"].iterkeys():
        md5s.append(md5)
    return md5s

def tmp_name():
    fd, temp_path = tempfile.mkstemp()
    os.close(fd)
    os.remove(temp_path)
    return temp_path

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
    parser.add_argument('--bigwigs', nargs='*')
    parser.add_argument('--labels', nargs='*')
    parser.add_argument('--md5s')
    parser.add_argument('--include')
    parser.add_argument('--exclude')
    parser.add_argument('--bin')
    parser.add_argument('--output')
    parser.add_argument('--assembly')
    args = parser.parse_args()
    if args.md5s == "None":
        args.md5s = []
    if args.bigwigs == ['None']:
        args.bigwigs = []
    if args.labels == ['None']:
        args.labels = []

    #read md5 list file for public data
    if args.md5s:
      md5_json = json.load(open(args.md5s))
      md5s = parse_md5s(md5_json)
    else:
      md5_json = {}
      md5s = []

    #public data paths
    include_path = config.REGION[args.assembly][args.include]
    exclude_path = config.REGION[args.assembly][args.exclude]

    #create temporary input files for geec executables
    input_list = []

    user_input_list = []
    for bw, label in itertools.izip(args.bigwigs, args.labels):
        user_hdf5 = tmp_name()
        user_filtered_hdf5 = tmp_name()
        user_input_list.append((bw, label, user_hdf5, user_filtered_hdf5))
        input_list.append((user_filtered_hdf5, label))


    public_path_dict = {}
    public_hdf5_paths_file = config.HDF5[args.assembly][args.bin][args.include][args.exclude]
    with open(public_hdf5_paths_file) as public_list:
        for line in public_list:
            line = line.split()
            public_path_dict[line[1]] = line[0]

    for md5 in md5s:
        if public_path_dict.get(md5, False):
            input_list.append((public_path_dict.get(md5), md5))

    correlation_file = tmp_name()
    input_list_path = create_input_list(input_list)

    # convert user bigwigs to hdf5 and filter it
    for raw_file, name, user_hdf5, user_filtered_hdf5 in user_input_list:
        to_hdf5(raw_file, name, args.assembly, user_hdf5, args.bin)
        filter_hdf5(name, args.assembly, user_hdf5, user_filtered_hdf5, args.bin, include_path, exclude_path)

    #correlate all uncorrelated matrix cells
    correlate(input_list_path, args.assembly, correlation_file, args.bin)

    #generate the final matrix
    make_matrix(input_list_path, args.assembly, correlation_file, args.output, args.md5s)

if __name__ == '__main__':
    main()
