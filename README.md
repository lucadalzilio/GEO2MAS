# stm 2d

Luca Dal Zilio, 2022, ETH Zurich

Reference:
Dal Zilio, L., van Dinther, Y., Gerya, T. V., and Pranger, C. C. (2018). 
Seismic behaviour of mountain belts controlled by plate convergence rate. 
Earth and Planetary Science Letters, 482, 81-92.

=====================================================

To run stm-2d on ETH Euler:

1. Load modules
	module load legacy
	module load intel
	module load hdf5

2. Export environmental variables
	export OMP_NUM_THREADS=6

3. Compile
	export I2STM_LIBS="-L/opt/intel/Compiler/11.1/075/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lifcore -mcmodel=medium -openmp -lhdf5"
	
	icc -std=c99 in2stm.c -o in2stm ${I2STM_LIBS} 
	
	icc -std=c99 i2stm.c -o i2stm ${I2STM_LIBS} 
	
4. Initialize the model setup
	./in2stm

5. Run 
	bsub -W 24:00 -n 6 -R 'rusage[mem=6144]' -J stm_t1 ./i2stm
	
	#if need to run longer than W for large-scale model use restarts
	
	bsub -W 24:00 -n 6 -R 'rusage[mem=6144]' -J stm_t1 -w "ended(stm_t1)" ./i2stm
