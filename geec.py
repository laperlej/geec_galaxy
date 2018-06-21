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
import multiprocessing

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), 'epigeec/epigeec/python/core'))
import main as epimain
import launcher
import make_matrix

PUBLIC_DATA_ROOT = "/geec-data/public"

import os.path

#directories
RESOURCE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "epigeec", "epigeec", "resource")
WIG_TO_BW = os.path.join(os.path.dirname(os.path.realpath(__file__)), "bin", "wigToBigWig")
MODULE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "geec_tools")

def analysis_path(script_name):
    return os.path.join(MODULE_DIR, 'geec_analysis', script_name)

GEEC_ANNOTATE = analysis_path('geec_annotate.py')
GEEC_ARI = analysis_path('geec_ari.py')
GEEC_SLICE = analysis_path('geec_slice.py')
GEEC_SLICE_FILE_NAME = analysis_path('geec_slice_file_name.py')

#chrom sizes
def chrom_sizes_path_maker(filename):
    return os.path.join(RESOURCE_DIR,'chrom_sizes',filename)

def get_chrom_sizes(assembly):
    assembly = assembly.lower()
    if assembly != 'saccer3':
        assembly = assembly + '.noy'
    else:
        assembly = assembly + '.can'
    filename = '{0}.chrom.sizes'.format(assembly)
    return chrom_sizes_path_maker(filename)

#regions
def filter_path_maker(filename):
    return os.path.join(RESOURCE_DIR, 'filter', filename)

def get_filter(assembly, content):
    if content == 'none':
        return filter_path_maker('none.bed')
    filename = "{0}.{1}.bed".format(assembly.lower(), content.lower())
    return filter_path_maker(filename)

#precalculated
def hdf5_path_maker(path):
    return os.path.join(PUBLIC_DATA_ROOT, path[0], path[1], path[2])

def get_resolution(num):
    to_human = {1:"1bp",
                10:"10bp",
                100: "100bp",
                1000: "1kb",
                10000: "10kb",
                100000: "100kb",
                1000000: "1mb",
                10000000: "10mb",
                100000000: "100mb"}
    return to_human[int(num)]

def get_matrix(assembly, resolution, include, exclude, metric="pearson"):
    filename = "{0}_{1}_{2}_{3}.mat".format(get_resolution(resolution), include, exclude, metric)
    path = [PUBLIC_DATA_ROOT, assembly, filename]
    return os.path.join(*path)

def get_hdf5(md5, assembly, resolution, include, exclude, metric="pearson"):
    ext = {"pearson":"hdf5",
           "spearman":"rank"}
    folder = "{0}_{1}_{2}".format(get_resolution(resolution), include, exclude)
    path = [assembly, folder, "{0}_{1}.{2}".format(md5, folder,ext[metric])]
    return hdf5_path_maker(path)

def to_hdf5(params):
    args, datatype, raw_file, name, user_hdf5, user_filtered_hdf5, include_path, exclude_path = params
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
        #continue
        return
    filter_hdf5(name, args.assembly, user_hdf5, user_filtered_hdf5, include_path, exclude_path)
    if args.metric == "spearman":
        rank_hdf5(user_filtered_hdf5)

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
    args = ["to_hdf5", "-bw", raw_file, get_chrom_sizes(assembly), resolution, user_hdf5]
    epimain.main(args)

def bg_to_hdf5(raw_file, name, assembly, user_hdf5, resolution):
    args = ["to_hdf5", "-bg", raw_file, get_chrom_sizes(assembly), resolution, user_hdf5]
    epimain.main(args)

def wig_to_bigwig(wig_file, bigwig_file):
    # wigToBigWig in.wig chrom.sizes out.bw
    wig = Wig(wig_file)
    wig.read()
    chromsizes_file = tmp_name()
    with open(chromsizes_file, "w") as chrom_size:
        chrom_size.write(str(wig))
    arguments = [WIG_TO_BW,
                 wig_file,
                 chromsizes_file,
                 bigwig_file]
    subprocess.call(arguments)


def filter_hdf5(name, assembly, user_hdf5, filtered_hdf5, include, exclude):
    args = ["filter"]
    if "all" not in include:
        args += ["--select", include]
    if "none" not in exclude:
        args += ["--exclude", exclude]
    args += [user_hdf5, get_chrom_sizes(assembly), filtered_hdf5]
    epimain.main(args)

def is_precalc(md5s, files, metric):
    # verify if nm
    return bool(md5s and metric == "pearson")

def correlate(input_list, assembly, mat_file):
    args = ["correlate", input_list, get_chrom_sizes(assembly), mat_file]
    epimain.main(args)

def correlate_nm(input_list1, input_list2, assembly, mat_file):
    print(False, input_list1, input_list2, get_chrom_sizes(assembly), mat_file)
    launcher.corr_nm(False, input_list1, input_list2, get_chrom_sizes(assembly), mat_file)

def launch_make_matrix(nn_mat_file, output_matrix, meta_json = ""):
    args = [nn_mat_file, output_matrix]
    if meta_json:
        args += ["-meta", meta_json]
    make_matrix.main(args)

def make_matrix_nm(nn_mat_file, nm_mat_file, precalc_matrix, output_matrix, meta_json = ""):
    args = ["-nm", nm_mat_file, "-mm", precalc_matrix, nn_mat_file, output_matrix]
    if meta_json:
        args += ["-meta", meta_json]
    make_matrix.main(args)

def slice_matrix(md5s, assembly, resolution, include, exclude, output):
    """
    python geec_slice_file_name.py matrix.mat fn_1 fn_2 fn_N > output.mat
    """
    arguments = ['python',
                 GEEC_SLICE_FILE_NAME,
                 get_matrix(assembly, resolution, include, exclude)]
    arguments += md5s.encode('ascii', 'ignore')
    print(arguments)
    with open(output, 'w') as output_file:
        subprocess.call(arguments, stdout=output_file)

def rank_hdf5(hdf5_path):
  h5f = h5py.File(hdf5_path, 'r+')
  for group in h5f:
    md5 = group
  for dset in h5f[md5]:
    data = h5f[md5][dset]
    n = len(data)
    assert(n < 3024616) #overflow on int64
    data[...] = scipy.stats.rankdata(data, method="ordinal")
    data.attrs['sumX'] = n * (n + 1) / 2 #sum of natural numbers
    data.attrs['sumXX'] = n * (n + 1) * (2 * n + 1) / 6 #sum of perfect squares
  h5f.close()

def create_input_list(input_list):
    file_path = tmp_name()
    with open(file_path, 'w') as input_file:
        for path, label in input_list:
            input_file.write('{0}\n'.format(path.strip()))
            #input_file.write('{0}\t{1}\n'.format(path.strip(),label.strip()))
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

    #public data paths
    include_path = get_filter(args.assembly, args.include)
    exclude_path = get_filter(args.assembly, args.exclude)
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
        hdf5_path = get_hdf5(md5, args.assembly, args.bin, args.include, args.exclude, args.metric)
        if os.path.isfile(hdf5_path):
            input_list2.append((hdf5_path, md5))
        else:
            print "{0} is missing".format(md5_json["datasets"][md5].get("file_name", "unknown"))

    mat_file_nn = tmp_name()
    mat_file_nm = tmp_name()
    mat_file_mm = tmp_name()

    # convert user bigwigs to hdf5 and filter it
    if user_input_list:
        p = multiprocessing.Pool(1)
        p_args = []
        for raw_file, datatype, name, user_hdf5, user_filtered_hdf5 in user_input_list:
            p_args.append((args, datatype, raw_file, name, user_hdf5, user_filtered_hdf5, include_path, exclude_path))
        p.map(to_hdf5, p_args)

    if is_precalc(md5s, args.files, args.metric):
        input_list_path1 = create_input_list(input_list1)
        input_list_path2 = create_input_list(input_list2)
        #correlate all uncorrelated matrix cells
        if input_list1:
            input_list_path1 = create_input_list(input_list1)
            input_list_path2 = create_input_list(input_list2)
            correlate(input_list_path1, args.assembly, mat_file_nn)
            correlate_nm(input_list_path1, input_list_path2, args.assembly, mat_file_nm)
            #generate the final matrix
            slice_matrix([x[1] for x in input_list2], args.assembly, args.bin, args.include, args.exclude, mat_file_mm)
            make_matrix_nm(mat_file_nn, mat_file_nm, mat_file_mm, args.output, args.md5s)
        else:
            input_list_path2 = create_input_list(input_list2)
            slice_matrix([x[1] for x in input_list2], args.assembly, args.bin, args.include, args.exclude, mat_file_nn)
            launch_make_matrix(mat_file_nn, args.output, args.md5s)
    else:
        input_list_path = create_input_list(input_list1 + input_list2)
        #correlate all uncorrelated matrix cells
        correlate(input_list_path, args.assembly, mat_file_nn)
        #generate the final matrix
        launch_make_matrix(mat_file_nn, args.output, args.md5s)
    #matrix_content = open(mat_file_nn).read()
    #for oldname, newname in :
    #with open(mat_file_nn, 'w') as f:
    #    f.write(matrix_content)


if __name__ == '__main__':
    main()
