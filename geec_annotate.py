import sys
import subprocess
import os.path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'geec_tools/python/utils'))
import config

def geec_annotate(args):
    arguments = [config.GEEC_ANNOTATE] + args

    subprocess.call(arguments)

def main():
    args = sys.argv[1:]
    geec_annotate(args)

if __name__ == '__main__':
    main()
