import sys
import subprocess
import os.path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'geec_tools/python/utils'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'geec_tools/geec_ari'))
import config
import geec_ari

def run_geec_ari(args):
    arguments = [config.GEEC_ARI] + args

    geec_ari.main(arguments)

def main():
    args = sys.argv[1:]
    run_geec_ari(args)

if __name__ == '__main__':
    main()
