PROG    := hello-world
PDA     := ../../pda_lib
LIB_OPT := pda
LIB_DEV := pda-dev

default: $(PROG)
all:     $(PROG) $(PROG)-dev

CXX = g++

CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic \
	-Winit-self -Wold-style-cast -Wcast-qual -Wundef \
	-Woverloaded-virtual \
	-Wconversion \
	-Wlogical-op \
	-Wzero-as-null-pointer-constant 

SANITIZE := -fsanitize=address \
	-fsanitize=undefined \
	-fsanitize=leak \
	-fsanitize=float-cast-overflow \
	-fno-omit-frame-pointer
	
CXXFLAGS_OPT := $(CXXFLAGS) -D NDEBUG -O3
CXXFLAGS_DEV := $(CXXFLAGS) $(SANITIZE) -O0 -g 

PDALIB_OPT  := $(PDA)/lib$(LIB_OPT).a
PDALIB_DEV  := $(PDA)/lib$(LIB_DEV).a

$(PROG): $(PROG).cpp $(PDALIB_OPT) makefile
	$(CXX) $(CXXFLAGS_OPT) -I$(PDA) -L$(PDA) $< -o $@ -l$(LIB_OPT)
$(PROG)-dev: $(PROG).cpp $(PDALIB_DEV) makefile
	$(CXX) $(CXXFLAGS_DEV) -I$(PDA) -L$(PDA) $< -o $@ -l$(LIB_DEV)

$(PDALIB_OPT):	
	cd $(dir $@); make $(notdir $@)
$(PDALIB_DEV):	
	cd $(dir $@); make $(notdir $@)

clean:
	rm -f $(PROG)
	rm -f $(PROG)-dev
	rm -f *.o
	rm -f *.a
	rm -f -r html
	rm -f *~
	rm -f *.stackdump
	rm -f *.exe
	rm -f core*
	cd $(PDA); make clean
