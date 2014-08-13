EXENAME		= run.x
SHELL		= /bin/bash
UNAME		= $(shell uname)
HOSTNAME	= $(shell hostname)
LOG		= compile.log

CC		= gcc
CFLAGS		= -g -O3 -fPIC

include Makefile.ifort
.SILENT:

LIBMODULESF90	= \
		fortpy.f90 \
		num_types.f90 \
		enumeration_types.f90 \
		numerical_utilities.f90 \
		vector_matrix_utilities.f90 \
		utilities.f90 \
		rational_mathematics.f90 \
		io_utils.f90 \
		symmetry_module.f90 \
		compare_structures.f90 \
		combinatorics.f90 \
		labeling_related.f90 \
		sorting.f90 \
		derivative_structure_generator.f90
MAINF90		= tester.f90
SRCF90		= $(LIBMODULESF90) $(MAINF90)
OBJSF90		= $(SRCF90:.f90=.o)

SRCC		= timing.c
HEADC		= timing.h
OBJSC		= timing.o

# Error handling
NEWFILE		= \#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#
ERR		= ******************************* ERROR *******************************
SHOW_LOG	= ( perl -pi -e 's/ [Ee]rror \#/\n\n\n$(ERR)\n*** error \#/' $(LOG); perl -pi -e 's/^\# 1 "/\n\n$(NEWFILE)\n\n\n/' $(LOG); grep -n -A3 -E "$(ERR)|$(NEWFILE)" $(LOG) )


all:	info $(EXENAME)

info: 
	echo -e "\nCompile time:" > $(LOG)
	date >> $(LOG)
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo "                     FORTPY"                           | tee -a $(LOG)
	echo "               >>> version 1.0 <<<                    "| tee -a $(LOG)         
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo -e "Compiling on system  : $(UNAME)"                    | tee -a $(LOG)
	echo -e "             machine : $(HOSTNAME)"                 | tee -a $(LOG)
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo -e "DEBUG mode\t:\t$(DEBUG)"                            | tee -a $(LOG)
	echo -e "MPI mode\t:\t$(UNCLE_MPI)"                          | tee -a $(LOG)
	echo -e "openMP mode\t:\t$(UNCLE_OMP)"                       | tee -a $(LOG)
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo "F90    : $(F90)"                                       | tee -a $(LOG)
	echo "FFLAGS : $(FFLAGS)"                                    | tee -a $(LOG)
	echo "CC     : $(CC)"                                        | tee -a $(LOG)
	echo "CFLAGS : $(CFLAGS)"                                    | tee -a $(LOG)
	echo "LDFLAGS: $(LDFLAGS)"                                   | tee -a $(LOG)
	echo "MKLpath:$(MKL)"                                        | tee -a $(LOG)
	echo "------------------------------------------------------"| tee -a $(LOG)
	echo ""                                                      | tee -a $(LOG)



$(EXENAME): $(OBJSF90) $(OBJSC)
	-rm $(EXENAME) 2> /dev/null
	echo -n "Linking... "
	-$(F90) $(LDFLAGS) -o $(EXENAME) $(OBJSC) $(OBJSF90) >> $(LOG) 2>> $(LOG)
	echo "done."
	if test -e $(EXENAME); then echo "Produced executable: $(EXENAME)"; else $(SHOW_LOG); echo "Error."; fi

$(OBJSF90): %.o: %.f90
	echo -n "Compiling: $^... "
	-$(F90) -c $(FFLAGS) $^ >> $(LOG) 2>> $(LOG)
	echo "done."

$(OBJSC): %.o: %.c $(HEADC)
	echo -n "Compiling: $^... "
	$(CC) $(CFLAGS) -c $^
	echo "done."
