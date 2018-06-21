export OMP_NUM_THREADS=10
export TCL_LIBRARY=/usr/share/tcl8.5
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/opt.usherbrooke.ca/CentOS6/szip/2.1/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/opt.usherbrooke.ca/intel/composerxe-2011.5.220/mkl/lib/intel64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/opt.usherbrooke.ca/intel/composerxe-2011.5.220/compiler/lib/intel64
source  /cvmfs/opt.usherbrooke.ca/lmod/lmod/init/profile
module use /cvmfs/opt.usherbrooke.ca/CentOS6/Modules/modulefiles.x86_64/
module use /cvmfs/opt.usherbrooke.ca/modulesfiles
module use /opt/Modules/modulefiles
module load python64/2.7.5 gcc/5.2.0 intel64/12.0.5.220 openmpi_intel64/1.4.3_ofed hdf5/1.8.16 boost64/1.61.0 cmake/2.8.8 R/3.0.1
