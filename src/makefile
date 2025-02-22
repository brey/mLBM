OBJS=lkcmod.o mlbm.o

#FLAGS=-fallow-argument-mismatch

mlbm:	$(OBJS)
	mpif90 -o mlbm $(OBJS) $(FLAGS)

lkcmod.o: lkcmod.f90
	mpif90 -c lkcmod.f90 $(FLAGS)
mlbm.o: mlbm.f90
	mpif90 -c mlbm.f90 $(FLAGS)
