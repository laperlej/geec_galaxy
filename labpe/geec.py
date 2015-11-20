import sys
import argparse
import subprocess
import tempfile
import itertools
import os

sys.path.insert(0, '/Users/Jon/Projects/geec_tools/python')
import config

def to_hdf5(input_list, assembly, user_hdf5, resolution):
    """Usage: to_hdf5 {input_list.txt}
                      {chrom_sizes}
                      {output.hdf5}
                      {bin_size}\n"""
    subprocess.call([config.TO_HDF5,
                     input_list,
                     config.CHROM_SIZE[assembly],
                     user_hdf5,
                     resolution])

def to_zscore(input_list, assembly, user_hdf5, user_zscore, resolution, include, exclude):
    """Usage: to_zscore {input_list.txt}
                        {chrom_sizes}
                        {input.hdf5}
                        {output.hdf5}
                        {bin_size}
                        {include.bed}
                        {exclude.bed}\n");"""
    subprocess.call([config.TO_ZSCORE,
                     input_list,
                     config.CHROM_SIZE[assembly],
                     user_hdf5,
                     user_zscore,
                     resolution,
                     include,
                     exclude
                    ])

def correlate(user_input_list, public_input_list, assembly, user_zscore, public_zscore, correlation_file, resolution):
    """Usage: correlation {input_list1}
                          {input_list2}
                          {chrom_sizes}
                          {input1.hdf5}
                          {input2.hdf5}
                          {output.results}
                          {bin_size}\n");"""
    subprocess.call([config.CORRELATION,
                     user_input_list,
                     public_input_list,
                     config.CHROM_SIZE[assembly],
                     user_zscore,
                     public_zscore,
                     correlation_file,
                     resolution
                     ])

def complete_matrix(user_input_list, public_input_list, assembly, correlation_file, public_matrix, output_matrix):
    """
    complete_matrix.py 
    chrom_sizes = sys.argv[1]
    list1 = sys.argv[2]
    list2 = sys.argv[3]
    correlation_file = sys.argv[4]
    matrix_file = sys.argv[5]
    output_matrix_path = sys.argv[6]
    """
    subprocess.call(['python', 
                     config.COMPLETE_MATRIX,
                     config.CHROM_SIZE[assembly],
                     user_input_list,
                     public_input_list,
                     correlation_file,
                     public_matrix,
                     output_matrix
                     ])

def create_input_list(paths, labels):
    file_path = tmp_name()
    with open(file_path, 'w') as input_file:
        for path, label in itertools.izip(paths, labels):
            input_file.write('{0}\t{1}\n'.format(path.strip(),label.strip()))
    return file_path

def parse_md5s(md5s_path):
    md5s = []
    if md5s_path:
        with open(md5s_path) as md5s_file:
            md5s_file.readline()
            for line in md5s_file:
                md5s.append(line)
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
    md5s = parse_md5s(args.md5s)

    #create temporary input files for geec executables
    user_input_list = create_input_list(args.bigwigs, args.labels)
    public_input_list = create_input_list(md5s, md5s)
    user_hdf5 = tmp_name()
    user_zscore = tmp_name()
    correlation_file = tmp_name()

    #public data paths
    public_zscore = config.ZSCORE[args.assembly][args.bin][args.include][args.exclude]
    public_matrix = config.MATRIX[args.assembly][args.bin][args.include][args.exclude]
    include_path = config.REGION[args.assembly][args.include]
    exclude_path = config.REGION[args.assembly][args.exclude]

    # convert user bigwigs to hdf5 and than zscore
    to_hdf5(user_input_list, args.assembly, user_hdf5, args.bin)
    to_zscore(user_input_list, args.assembly, user_hdf5, user_zscore, args.bin, include_path, exclude_path)

    #correlate all uncorrelated matrix cells
    correlate(user_input_list, public_input_list, args.assembly, user_zscore, public_zscore, correlation_file, args.bin)

    #generate the final matrix
    complete_matrix(user_input_list, public_input_list, args.assembly, correlation_file, public_matrix, args.output)

if __name__ == '__main__':
    main()
