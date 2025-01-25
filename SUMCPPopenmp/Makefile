
exe=sum.x

all: ${exe}

.PHONY: all clean allclean

##################################
# setup for CfCA XD2000 intel env begin
##################################

cc= mpiicpx
coptopenmp = -fiopenmp
copt = -O2

# setup for CfCA XD2000 end
##################################

.SUFFIXES: .cpp .o

.cpp.o: 
	${cc} ${copt} ${coptopenmp} -c $< 

OBJS=main.o mpiroutines.o

mpiroutines.o:
main.o: main.cpp mpiroutines.o 

${exe}: ${OBJS}
	${cc} ${copt} ${coptopenmp} ${OBJS} -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *~
