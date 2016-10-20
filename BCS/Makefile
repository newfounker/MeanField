include make.inc

OBJS= 

all: version compile completion
debug: debug version compile completion

debug: FFLAG=$(DFLAG)

compile: $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FFLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)$(BRANCH)

.f90.o:	
	$(FC) $(FFLAG) -c $< 

completion:
	scifor_completion.sh $(DIR)/$(EXE).f90

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)

