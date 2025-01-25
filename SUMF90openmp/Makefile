
exe=sum.x

all: ${exe}

.PHONY: all clean allclean

##################################
# setup for CfCA XD2000 intel env begin
##################################

fc= mpiifx -fc=ifort -diag-disable=10448
foptopenmp = -qopenmp
fopt = -g -traceback -O2

# setup for CfCA XD2000 end
##################################

##################################
# setup for CfCA XD2000 intel(ifx) env begin
##################################

#fc= mpiifx
#foptopenmp = -fiopenmp
#fopt = -g -traceback -O2 -qopt-report=2
#foptlink= -Wl,--no-relax
#fopt=-g -traceback -check all -fpe0

# setup for CfCA XD2000 end
##################################

##################################
# setup for Fugaku
##################################
#fc= mpifrtpx -Free
#foptopenmp = -Kopenmp
#fopt = -Kfast
# setup for Fugaku end
##################################


.SUFFIXES: .f90 .o

.f90.o: 
	${fc} ${fopt} ${foptopenmp} -c $< 

OBJS=main.o mpiroutines.o

mpiroutines.o:
main.o: main.f90 mpiroutines.o 

${exe}: ${OBJS}
	${fc} ${fopt} ${foptopenmp} ${OBJS} ${foptlink} -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~
