import sys
import subprocess
import os.path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'geec_tools/python/utils'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'geec_tools/geec_analysis'))
import config
import geec_annotate

def run_geec_annotate(args):
    arguments = [config.GEEC_ANNOTATE] + args

    geec_annotate.main(arguments)

def main():
    args = sys.argv[1:]
    run_geec_annotate(args)

if __name__ == '__main__':
    main()
