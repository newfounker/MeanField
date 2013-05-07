include sfmake.inc
FC=ifort
HERE   = `pwd`
EXE=ahmbcs_bethe
DIREXE   = $(HERE)/test

OBJS     = 
OBJS_DEB =

code: 	$(OBJS)
	@echo "compile: standard"
	${FC} ${STD} ${OBJS} ${EXE}.f90 -o ${DIREXE}/${EXE} ${SFMODS} ${SFLIBS} ${MORELIB}
	@echo ""

opt: 	$(OBJS)
	@echo "compile: optimized"
	${FC} ${OPT}  ${OBJS} ${EXE}.f90 -o ${DIREXE}/${EXE} ${SFMODS} ${SFLIBS}
	@echo ""

debug:	$(OBJS_DEB)	
	@echo "compile: debug"
	${FC} ${DEB} ${OBJS} ${EXE}.f90 -o ${DIREXE}/${EXE}  ${SFMODS} ${SFLIBS} ${MORELIB}
	@echo ""

clean: 
	@echo "Cleaning:"
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~
	@rm -f ${DIREXE}/${EXE}
	@echo ""
