CC1 = mpicc
CC2 = gcc
FLAGS = -lm


all:
	$(CC1) *.c -D MIGSELECTRELEASE -o MigSelect $(FLAGS)

singlecpu:
	$(CC2) *.c -U MPI_ENABLED -D MIGSELECTRELEASE -o MigSelect_singlecpu $(FLAGS)

clean:
	rm MigSelect MigSelect_singlecpu
	rm *.o
