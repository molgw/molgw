# This file is part of MOLGW
# Author: Fabien Bruneval

-include ../src/my_machine.arch
-include ../my_machine.arch
-include ./my_machine.arch

.PHONY: clean

basis_2_basisrel.x: basis_2_basisrel.cpp
	$(CXX) $(CXXFLAGS) -c basis_2_basisrel.cpp
	$(CXX) basis_2_basisrel.cpp -o ../basis_2_basisrel.x


clean:
	$(RM) ../basis_2_basisrel.x basis_2_basisrel.o
