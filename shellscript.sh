make clean
make
mpirun -np 64 ./main inputs/Silica_beta_cristobalite.dat output/data/densityandenergy.dat output/data/RDF.dat output/data/BAD.dat
make -C output/gnuplot/
