TARGET = a.out
OBJECTS = GaussLegen_module.o fct.o main.o
F90 = ifort
COMMON_MOD = GaussLegen_module.f90
RM := rm

${TARGET} : ${OBJECTS}
	${F90}  ${OBJECTS}

${OBJECTS} : ${COMMON_MOD} fct.f90 main.f90
	${F90} -c -CB fct.f90  GaussLegen_module.f90 main.f90

.PHONY : clean
clean :
	${RM} *.o *.mod *~
