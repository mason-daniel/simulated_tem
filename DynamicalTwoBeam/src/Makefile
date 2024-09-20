.SUFFIXES: .o .F .f .f90 .c .java .class .F90

#	all the source in dependency order


IRSIDE_DEF = \
inputREAL.o         \

ITSIDE_DEF = \
inputTDS.o         \

RFSIDE_DEF = \
REALfit.o         \

TFSIDE_DEF = \
TDSfit.o         \

DWSIDE_DEF = \
DebyeWallerFactor.o         \

CSIDE_DEF = \
NBAX_StringTokenizers.o         \
Lib_CommandLineArguments.o      \
Lib_RelativisticElectrons.o     \
Lib_CrystalStructureFactor.o     \
crystalStructureFactor.o         \

################################

include Makefile.inc

IRSIDE = inputREAL.exe

ITSIDE = inputTDS.exe

RFSIDE = REALfit.exe

TFSIDE = TDSfit.exe

DWSIDE = DebyeWallerFactor.exe

CSIDE = crystalStructureFactor.exe


all:
#	cp -u ${MYF90LIB}/NBAX_StringTokenizers.f90 .
#	cp -u ${MYF90LIB}/Lib_CommandLineArguments.f90 .
#	cp -u ${MYF90LIB}/Lib_RelativisticElectrons.f90 .	
	make irside
	make itside
	make rfside
	make tfside
	make dwside
	make cside

irside:  $(IRSIDE)

itside:  $(ITSIDE)

rfside:  $(RFSIDE)

tfside:  $(TFSIDE)

dwside:  $(DWSIDE)

cside:  $(CSIDE)



$(IRSIDE) : $(IRSIDE_DEF)
	$(LDR) $(LFLAGS) $(FFLAGS) -o $(IRSIDE) $(IRSIDE_DEF) $(LIBS)
	cp -u $(IRSIDE) ../bin/.

$(ITSIDE) : $(ITSIDE_DEF)
	$(LDR) $(LFLAGS) $(FFLAGS) -o $(ITSIDE) $(ITSIDE_DEF) $(LIBS)
	cp -u $(ITSIDE) ../bin/.

$(RFSIDE) : $(RFSIDE_DEF)
	$(LDR) $(LFLAGS) $(FFLAGS) -o $(RFSIDE) $(RFSIDE_DEF) $(LIBS)
	cp -u $(RFSIDE) ../bin/.

$(TFSIDE) : $(TFSIDE_DEF)
	$(LDR) $(LFLAGS) $(FFLAGS) -o $(TFSIDE) $(TFSIDE_DEF) $(LIBS)
	cp -u $(TFSIDE) ../bin/.

$(DWSIDE) : $(DWSIDE_DEF)
	$(LDR) $(LFLAGS) $(FFLAGS) -o $(DWSIDE) $(DWSIDE_DEF) $(LIBS)
	cp -u $(DWSIDE) ../bin/.

$(CSIDE) : $(CSIDE_DEF)
	$(LDR) $(LFLAGS) $(FFLAGS) -o $(CSIDE) $(CSIDE_DEF) $(LIBS)
	cp -u $(CSIDE) ../bin/.


#	rules for compiling source


.F90.o:
	$(FC90) -c $(FFLAGS) $*.F90

.f90.o:
	$(FC90) -c $(FFLAGS) $*.f90

.f.o:
	$(FCold) -c $(FFLAGS) $*.f	 

.F.o:
	$(FC) -c $(FFLAGS) $*.F -FR

.c.o:
	$(CC) -c $(CFLAGS) $*.c

.java.class:
	$(JAVAC) $(JFLAGS) $*.java

#	rules for tidying up


clean:
	rm -f *.o *.run *.mod *.il ifc* *.a *.so *.mpi *.class


#################
