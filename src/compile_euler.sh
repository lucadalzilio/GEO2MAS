# load modules
module load legacy
module load intel
module load hdf5

# export environmental variables
export OMP_NUM_THREADS=6

# compile
# if want to debug, add -g after icc
export I2STM_LIBS="-L/opt/intel/Compiler/11.1/075/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lifcore -mcmodel=medium -openmp -lhdf5"
icc -std=c99 in2stm.c -o in2stm ${I2STM_LIBS} 
icc -std=c99 i2stm.c -o i2stm ${I2STM_LIBS} 

./in2stm

