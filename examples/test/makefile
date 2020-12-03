PDA     := ../../pda_lib
LIBNAME := pda_deb
#LIBNAME := pda
PROG    := test

CXX = g++

SANITIZE := -fsanitize=address \
	-fsanitize=undefined \
	-fsanitize=leak \
	-fsanitize=float-cast-overflow \
	-fno-omit-frame-pointer

CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -g -O0 \
	-Winit-self -Wold-style-cast -Wcast-qual -Wundef \
	-Woverloaded-virtual \
	-Wconversion \
	-Wlogical-op \
	-Wzero-as-null-pointer-constant \
	$(SANITIZE)

PDALIB	=$(PDA)/lib$(LIBNAME).a

$(PROG): $(PROG).cpp $(PDALIB) makefile
	$(CXX) $(CXXFLAGS) -I$(PDA) -L$(PDA) $< -o $@ -l$(LIBNAME)

$(PDALIB):	
	cd $(dir $(PDALIB)); make $(notdir $(PDALIB))

clean:	makefile
	rm -f $(PROG)
	rm -f *.o
	rm -f *.a
	rm -f -r html
	rm -f *~
	rm -f *.stackdump
	rm -f *.exe
	rm -f core*
	cd $(dir $(PDALIB)); make clean
