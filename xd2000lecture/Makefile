
exe_mpi=kh_mpi.x
exe_hybrid=kh_hybrid.x
exe_openmp=kh_openmp.x

all: mpi

.PHONY: all mpi hybrid openmp clean

mpi: ${exe_mpi}

hybrid: ${exe_hybrid}

openmp: ${exe_openmp}

##################################
# setup for CfCA XD2000 cpe env (gnu) begin
##################################
fc= ftn
fopt = -O2

# setup for CfCA XD2000 end
##################################

##################################
# setup for CfCA XD2000 cpe env (cray) begin
##################################

#fc= ftn 
#fopt = -O2

# setup for CfCA XD2000 end
##################################

##################################
# setup for CfCA XD2000 cpe env (intel) begin
##################################
#fc= ftn 
#fopt = -O2

# setup for CfCA XD2000 end
##################################


##################################
# setup for CfCA XD2000 intel env begin
##################################

#fc= mpiifx -fc=ifort -diag-disable=10448
#foptopenmp = -qopenmp
#fopt = -g -traceback -O2

# setup for CfCA XD2000 end
##################################

##################################
# setup for CfCA XD2000 intel(ifx) env begin
##################################

#fc= mpiifx
#foptopenmp = -fiopenmp
#fopt = -g -traceback -O2 # -qopt-report=2
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


${exe_mpi}: main_mpi.f90 mpiroutines.f90
	${fc} ${fopt} mpiroutines.f90 main_mpi.f90 -o ${exe_mpi}

${exe_hybrid}: main_mpi.f90 mpiroutines.f90
	${fc} ${fopt} -fopenmp mpiroutines.f90 main_mpi.f90 -o ${exe_hybrid}

${exe_openmp}: main_openmp.f90 mpiroutines.f90
	${fc} ${fopt} -fopenmp mpiroutines.f90 main_openmp.f90 -o ${exe_openmp}

#################
# clean up
#################
clean:
	rm -f ${exe_mpi} ${exe_hybrid} ${exe_openmp} *.o *.mod *~
