#########################################################

# REPLACE the main_prg_name below with the name of
# your main program without the ".f90" extension
# REPLACE the module_name below with your module
# name without the ".f90" extension

# PROGRAM NAME (no .f90)
main=leapfrog

# MODULE NAME  (no .f90)
mod1=shared_data
mod2=diagnostics
mod3=control
mod4=setup
mod5=update_variables

#########################################################

.PHONY: clean

cmplr=gfortran

# flags=-O0 -g -Wall -Wextra -pedantic -fbounds-check \
#       -ffpe-trap=invalid,zero,overflow -Wno-unused-parameter \
# 			-Wrealloc-lhs-all -fno-realloc-lhs
flags = -O3

objects = $(mod1).o $(mod2).o $(mod3).o $(mod4).o $(mod5).o $(main).o

$(main).exe   : $(objects)
	$(cmplr) $(flags) -o $(main).exe $(objects)

$(mod1).o    : $(mod1).f90
	$(cmplr) $(flags) -c $(mod1).f90

$(mod2).o    : $(mod2).f90
	$(cmplr) $(flags) -c $(mod2).f90

$(mod3).o    : $(mod3).f90
	$(cmplr) $(flags) -c $(mod3).f90

$(mod4).o    : $(mod4).f90
	$(cmplr) $(flags) -c $(mod4).f90

$(mod5).o    : $(mod5).f90
	$(cmplr) $(flags) -c $(mod5).f90

$(main).o  : $(main).f90 $(mod1).f90 $(mod2).f90 $(mod3).f90 $(mod4).f90 $(mod5).f90
	$(cmplr) $(flags) -c $(main).f90

clean:
	-rm -f *.o *.mod *.exe
	-rm -f Data/*
	-rm -rf Figures/*
