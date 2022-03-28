# load modules
module load legacy
module load intel
module load hdf5

# export environmental variables
export OMP_NUM_THREADS=6

# run 
#bsub -I -W 00:10 -n 8 -R 'rusage[mem=6144]' -J stm_t1 ./i2stm

# run 
bsub -W 16:00 -n 6 -R 'rusage[mem=6144]' -J stm_t1 ./i2stm
# if need to run longer than W for large-scale model use restarts
bsub -W 16:00 -n 6 -R 'rusage[mem=6144]' -J stm_t1 -w "ended(stm_t1)" ./i2stm
bsub -W 16:00 -n 6 -R 'rusage[mem=6144]' -J stm_t1 -w "ended(stm_t1)" ./i2stm
#bsub -W 36:00 -n 6 -R 'rusage[mem=6144]' -J stm_t1 -w "ended(stm_t1)" ./i2stm
#bsub -W 36:00 -n 6 -R 'rusage[mem=6144]' -J stm_t1 -w "ended(stm_t1)" ./i2stm
#bsub -W 36:00 -n 6 -R 'rusage[mem=6144]' -J stm_t1 -w "ended(stm_t1)" ./i2stm
#bsub -W 36:00 -n 6 -R 'rusage[mem=6144]' -J stm_t1 -w "ended(stm_t1)" ./i2stm
