LOG		= compile.log
#FFLAGS		= -c
.SILENT:

# sources
WSRCF90		=  \
		f90wrap_fortpy.f90 \

SRCF90		=  \
		fortpy.f90 \

OBJSF90		= $(SRCF90:.f90=.o)

all: $(OBJSF90) $(WSRCF90)
	-rm _ifortpy.so
	echo "Creating shared library _ifortpy.so"
	-f2py-f90wrap -c -m _ifortpy -I. $(OBJSF90) $(WSRCF90) -I. --f90flags="-fmax-stack-var-size=1024">> f2py.log
	echo "done."

$(OBJSF90): %.o: %.f90
	echo "Compiling: $^"
	-$(F90) -c $^
	-$(F90) -c $(FFLAGS) $^ -o $@ >> $(LOG) 2>> $(LOG)
	echo "done."

$(WSRCF90): 
	f90wrap -k kind_map -m ifortpy $(SRCF90)

clean:
	-rm *.o *.mod *.so f90wrap* f2py.log
