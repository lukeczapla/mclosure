# this is the Makefile, type "make" to install the programs in mclosure/bin
# code written by Luke Czapla (Rutgers University)
#

CC=/usr/bin/g++

INSTDIR=../bin

CFLAGS=-O -Wall -D_FILE_OFFSET_BITS=64

LIBS=m

all: closuremc-lac closuremc closuremc-proteins closuremc-HU-hinge out2pdb closuremc-hinge radial

clean: 
	rm *.o

closuremc: closuremc.o matrix.o dna.o translate_tools.o bdna.o closure.o spatial.o histogram.o file.o trajectory.o proteinPDB.o
	$(CC) $(CFLAGS) -o $(INSTDIR)/closuremc closuremc.o translate_tools.o matrix.o dna.o bdna.o closure.o spatial.o histogram.o file.o trajectory.o proteinPDB.o -l$(LIBS)

closuremc-lac: closuremc-lac.o matrix.o dna.o bdna.o closure-lac.o spatial.o histogram.o file.o trajectory.o proteinPDB.o translate_tools.o
	$(CC) $(CFLAGS) -o $(INSTDIR)/closuremc-lac closuremc-lac.o matrix.o dna.o bdna.o closure-lac.o spatial.o histogram.o file.o trajectory.o translate_tools.o proteinPDB.o -l$(LIBS)

radial: radial.o matrix.o dna.o bdna.o closure-radial.o spatial.o histogram.o file.o trajectory.o proteinPDB.o
	$(CC) $(CFLAGS) -o $(INSTDIR)/radial radial.o matrix.o dna.o bdna.o closure-radial.o spatial.o histogram.o file.o trajectory.o proteinPDB.o

closuremc-parallel: closuremc-parallel.o matrix.o dna.o bdna.o closure-parallel.o spatial.o histogram.o file.o trajectory.o translate_tools.o proteinPDB.o
	$(CC) $(CFLAGS) -lpthread -o $(INSTDIR)/closuremc-parallel closuremc-parallel.o matrix.o dna.o bdna.o closure-parallel.o spatial.o histogram.o file.o trajectory.o translate_tools.o proteinPDB.o -l$(LIBS)

closuremc-pparallel: closuremc-pparallel.o matrix.o dna.o bdna.o closure-pparallel.o spatial.o histogram.o file.o trajectory.o translate_tools.o proteinPDB.o
	$(CC) $(CFLAGS) -lpthread -o $(INSTDIR)/closuremc-pparallel closuremc-pparallel.o matrix.o dna.o bdna.o closure-pparallel.o spatial.o histogram.o file.o trajectory.o translate_tools.o proteinPDB.o -l$(LIBS)


closuremc-HU: closuremc-HU.o matrix.o dna.o bdna.o closure-HU.o spatial.o histogram.o file.o trajectory.o
	$(CC) $(CFLAGS) -o $(INSTDIR)/closuremc-HU closuremc-HU.o matrix.o dna.o bdna.o closure-HU.o spatial.o histogram.o file.o trajectory.o -l$(LIBS)

closuremc-proteins: closuremc-proteins.o matrix.o dna.o bdna.o closure-proteins.o spatial.o histogram.o file.o trajectory.o translate_tools.o proteinPDB.o
	$(CC) $(CFLAGS) -o $(INSTDIR)/closuremc-proteins closuremc-proteins.o matrix.o dna.o bdna.o closure-proteins.o spatial.o histogram.o file.o trajectory.o translate_tools.o proteinPDB.o -l$(LIBS)

closuremc-HU-end: closuremc-HU-end.o matrix.o dna.o bdna.o closure-HU-end.o spatial.o histogram.o file.o trajectory.o translate_tools.o
	$(CC) $(CFLAGS) -o $(INSTDIR)/closuremc-HU-end closuremc-HU-end.o matrix.o dna.o bdna.o closure-HU-end.o spatial.o histogram.o file.o trajectory.o translate_tools.o -l$(LIBS)

closuremc-HU-long: closuremc-HU-long.o matrix.o dna.o bdna.o closure-HU-long.o spatial.o histogram.o file.o trajectory.o translate_tools.o
	$(CC) $(CFLAGS) -o $(INSTDIR)/closuremc-HU-long closuremc-HU-long.o matrix.o dna.o bdna.o closure-HU-long.o spatial.o histogram.o file.o trajectory.o translate_tools.o -l$(LIBS)

closuremc-disk: closuremc-disk.o matrix.o dna.o bdna.o closure-disk.o spatial.o histogram.o file.o trajectory.o
	$(CC) $(CFLAGS) -o $(INSTDIR)/closuremc-disk closuremc-disk.o matrix.o dna.o bdna.o closure-disk.o spatial.o histogram.o file.o trajectory.o -l$(LIBS)

closuremc-hinge: closuremc-hinge.o matrix.o dna.o bdna.o closure-hinge.o spatial.o histogram.o file.o trajectory.o
	$(CC) $(CFLAGS) -o $(INSTDIR)/closuremc-hinge closuremc-hinge.o matrix.o dna.o bdna.o closure-hinge.o spatial.o histogram.o file.o trajectory.o -l$(LIBS)

closuremc-HU-hinge: closuremc-HU-hinge.o matrix.o dna.o bdna.o closure-HU-hinge.o spatial.o histogram.o file.o trajectory.o
	$(CC) $(CFLAGS) -o $(INSTDIR)/closuremc-HU-hinge closuremc-HU-hinge.o matrix.o dna.o bdna.o closure-HU-hinge.o spatial.o histogram.o file.o trajectory.o -l$(LIBS)


closuremc-HU-hinge.o: closuremc-HU-hinge.cpp
	$(CC) $(CFLAGS) -c closuremc-HU-hinge.cpp

out2pdb: out2pdb.o matrix.o dna.o bdna.o
	$(CC) -o $(INSTDIR)/out2pdb out2pdb.o matrix.o dna.o bdna.o

matrix.o: matrix/matrix.cpp
	$(CC) $(CFLAGS) -c matrix/matrix.cpp 

out2pdb.o: out2pdb.cpp
	$(CC) $(CFLAGS) -c out2pdb.cpp

file.o: file.cpp
	$(CC) $(CFLAGS) -c file.cpp

dna.o: dna/dna.cpp
	$(CC) $(CFLAGS) -c dna/dna.cpp

bdna.o: bdna/bdna.cpp
	$(CC) $(CFLAGS) -c bdna/bdna.cpp

proteinPDB.o: proteinPDB/proteinPDB.cpp
	$(CC) $(CFLAGS) -c proteinPDB/proteinPDB.cpp

closure.o: closure/closure.cpp
	$(CC) $(CFLAGS) -c closure/closure.cpp

closure-radial.o: closure/closure-radial.cpp
	$(CC) $(CFLAGS) -c closure/closure-radial.cpp

closure-parallel.o: closure/closure-parallel.cpp
	$(CC) $(CFLAGS) -c closure/closure-parallel.cpp

closure-pparallel.o: closure/closure-pparallel.cpp
	$(CC) $(CFLAGS) -c closure/closure-pparallel.cpp

closure-lac.o: closure/closure-lac.cpp
	$(CC) $(CFLAGS) -c closure/closure-lac.cpp

closure-HU.o: closure/closure-HU.cpp
	$(CC) $(CFLAGS) -c closure/closure-HU.cpp

closure-HU-long.o: closure/closure-HU-long.cpp
	$(CC) $(CFLAGS) -c closure/closure-HU-long.cpp

closure-HU-hinge.o: closure/closure-HU-hinge.cpp
	$(CC) $(CFLAGS) -c closure/closure-HU-hinge.cpp

closure-HU-end.o: closure/closure-HU-end.cpp
	$(CC) $(CFLAGS) -c closure/closure-HU-end.cpp

closure-proteins.o: closure/closure-proteins.cpp
	$(CC) $(CFLAGS) -c closure/closure-proteins.cpp

closure-disk.o: closure/closure-disk.cpp
	$(CC) $(CFLAGS) -c closure/closure-disk.cpp

closure-hinge.o: closure/closure-hinge.cpp
	$(CC) $(CFLAGS) -c closure/closure-hinge.cpp

closure-HU-hinge.o: closure/closure-HU-hinge.cpp
	$(CC) $(CFLAGS) -c closure/closure-HU-hinge.cpp

spatial.o: spatial/spatial.cpp
	$(CC) $(CFLAGS) -c spatial/spatial.cpp

energy.o: energy/energy.cpp
	$(CC) $(CFLAGS) -c energy/energy.cpp

plot.o: plot/plot.cpp
	$(CC) $(CFLAGS) -c plot/plot.cpp

translate_tools.o: translate_tools.cpp
	$(CC) $(CFLAGS) -c translate_tools.cpp

histogram.o: histogram/histogram.cpp
	$(CC) $(CFLAGS) -c histogram/histogram.cpp

radial.o: radial.cpp
	$(CC) $(CFLAGS) -c radial.cpp

trajectory.o: trajectory/trajectory.cpp
	$(CC) $(CFLAGS) -c trajectory/trajectory.cpp

closuremc-lac.o: closuremc-lac.cpp
	$(CC) $(CFLAGS) -c closuremc-lac.cpp

closuremc.o: closuremc.cpp
	$(CC) $(CFLAGS) -c closuremc.cpp

closuremc-parallel.o: closuremc-parallel.cpp
	$(CC) $(CFLAGS) -c closuremc-parallel.cpp

closuremc-pparallel.o: closuremc-pparallel.cpp
	$(CC) $(CFLAGS) -c closuremc-pparallel.cpp

closuremc-hinge.o: closuremc-hinge.cpp
	$(CC) $(CFLAGS) -c closuremc-hinge.cpp

