import sys
import subprocess
import os.path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'geec_tools/python/utils'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'geec_tools/geec_analysis'))
import config
import geec_slice

def run_geec_slice(args):
    arguments = [config.GEEC_SLICE] + args

    geec_slice.main(arguments)

def main():
    args = sys.argv[1:]
    run_geec_slice(args)

if __name__ == '__main__':
    main()
