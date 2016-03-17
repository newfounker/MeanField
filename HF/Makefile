#include sfmake.inc
FC=ifort
HERE   = `pwd`
EXE=hmhf_1d
DIREXE   = $(HOME)/.bin

OBJS     = 
OBJS_DEB =
LIBS=-ldmftt -lscifor -llapack -lblas -lquadpack -lminpack

code: 	$(OBJS)
	@echo "compile: standard"
	${FC} ${FFLAG} ${OBJS} ${EXE}.f90 -o ${DIREXE}/${EXE}  ${LIBS} 
	@echo ""


clean: 
	@echo "Cleaning:"
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~
	@rm -f ${DIREXE}/${EXE}
	@echo ""
