#!/bin/bash
cd ~/PAP_Lab4_MPI_LatticeBoltzmann
make
# mpirun -np 8 ./check_comm -e 4 -p rank
# mpirun -np 8 ./check_comm -e 4 -p position -s both
# mpirun -np 1 ./lbm -e 0
# mpirun -np 8 ./lbm -e 4
mpirun --mca pml ob1 --mca btl tcp,self -np 8 ./check_comm -e 4 -p rank
mpirun --mca pml ob1 --mca btl tcp,self -np 8 ./check_comm -e 4 -p position -s both
mpirun --mca pml ob1 --mca btl tcp,self -np 1 ./lbm -e 0
./display --gnuplot output.raw 3 | md5sum > ref.md5
mpirun --mca pml ob1 --mca btl tcp,self -np 8 ./lbm -e 4
./display --gnuplot output.raw 3 | md5sum -c ref.md5
./gen_animate_gif.sh output.raw output.gif
