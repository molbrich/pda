PDA     := ../../pda_lib
LIBNAME := pda
PROG    := rc

PDALIB  := $(PDA)/lib$(LIBNAME).a

$(PROG): $(PROG).cpp $(PDALIB) makefile
	g++ -pedantic -g -I$(PDA) -L$(PDA) $< -o $@ -l$(LIBNAME)

$(PDALIB):	
	cd $(dir $(PDALIB)); make $(notdir $(PDALIB))

clean:	makefile
	rm -f $PROG$
	rm -f *.o
	rm -f *.a
	rm -f -r html
	rm -f *~
	rm -f *.stackdump
	rm -f *.exe
	rm -f core*
	cd $(dir $(PDALIB)); make clean
