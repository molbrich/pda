SRC	= pda_pda.cpp pda_pdv.cpp pda_matrix.cpp pda_vector.cpp pda_util.cpp pda_powersiterator.cpp
INC	= pda_pda.h pda_pdv.h pda_matrix.h pda_vector.h pda_util.h pda_powersiterator.h

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

all:	libpda.a libpda_deb.a

pda_pda.o:	pda_pda.cpp pda_pda.h makefile
	$(CXX) $(CXXFLAGS) -g -c $< -o $@
pda_pdv.o:	pda_pdv.cpp pda_pdv.h makefile
	$(CXX) $(CXXFLAGS) -g -c $< -o $@
pda_matrix.o:	pda_matrix.cpp pda_pda.h pda_pdv.h pda_matrix.h pda_vector.h makefile
	$(CXX) $(CXXFLAGS) -g -c $< -o $@
pda_vector.o:	pda_vector.cpp pda_pda.h pda_pdv.h pda_matrix.h pda_vector.h makefile
	$(CXX) $(CXXFLAGS) -g -c $< -o $@
pda_util.o:	pda_util.cpp pda_pda.h pda_pdv.h pda_matrix.h pda_vector.h pda_util.h makefile
	$(CXX) $(CXXFLAGS) -g -c $< -o $@
pda_powersiterator.o:	pda_powersiterator.cpp pda_pda.h pda_pdv.h pda_matrix.h pda_vector.h pda_util.h pda_util.cpp makefile
	$(CXX) $(CXXFLAGS) -g -c $< -o $@

libpda_deb.a:	pda_pda.o pda_pdv.o pda_matrix.o pda_vector.o pda_util.o pda_powersiterator.o $(INC) makefile
	ar -rcv libpda_deb.a pda_pda.o pda_pdv.o pda_matrix.o pda_vector.o pda_util.o pda_powersiterator.o
	ranlib libpda_deb.a
	
CXXFLAGS_OPT = -pedantic -D NDEBUG -O3 

pda_pda_opt.o:	pda_pda.cpp pda_pda.h makefile
	$(CXX) $(CXXFLAGS_OPT) -g -c $< -o $@
pda_pdv_opt.o:	pda_pdv.cpp pda_pdv.h pda_pda.h makefile
	$(CXX) $(CXXFLAGS_OPT) -g -c $< -o $@
pda_matrix_opt.o:	pda_matrix.cpp pda_pda.h pda_pdv.h pda_matrix.h pda_vector.h makefile
	$(CXX) $(CXXFLAGS_OPT) -g -c $< -o $@
pda_vector_opt.o:	pda_vector.cpp pda_pda.h pda_pdv.h pda_matrix.h pda_vector.h makefile
	$(CXX) $(CXXFLAGS_OPT) -g -c $< -o $@
pda_util_opt.o:	pda_util.cpp pda_pda.h pda_pdv.h pda_matrix.h pda_vector.h pda_util.h makefile
	$(CXX) $(CXXFLAGS_OPT) -g -c $< -o $@
pda_powersiterator_opt.o:	pda_powersiterator.cpp pda_pda.h pda_pdv.h pda_matrix.h pda_vector.h pda_util.h pda_util.cpp makefile
	$(CXX) $(CXXFLAGS_OPT) -g -c $< -o $@

libpda.a:	pda_pda_opt.o pda_pdv_opt.o pda_matrix_opt.o pda_vector_opt.o pda_util_opt.o pda_powersiterator_opt.o $(INC) makefile
	ar -rcv libpda.a pda_pda_opt.o pda_pdv_opt.o pda_matrix_opt.o pda_vector_opt.o pda_util_opt.o pda_powersiterator_opt.o
	ranlib libpda.a

.PHONY: doc, clean
doc:	makefile $(SRC) $(INC)
	doxygen doxyConfig

clean:	makefile
	rm -rf doc
	rm -f main
	rm -f test
	rm -f test_m32
	rm -f testcheck
	rm -f test_check
	rm -f test_new
	rm -f test_new_m32
	rm -f example
	rm -f *.o
	rm -f *.a
	rm -f -r html
	rm -f *~
	rm -f *.stackdump
	rm -f *.exe
	rm -f core*
