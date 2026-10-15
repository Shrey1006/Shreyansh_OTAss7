CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wno-sign-compare -Wno-unused-variable
all: moo
moo: main.cpp moo_utils.cpp moo_utils.hpp problems.hpp classical.hpp evo.hpp swarm.hpp aco.hpp grad.hpp
	$(CXX) $(CXXFLAGS) main.cpp moo_utils.cpp -o moo
test: moo
	./moo unit
clean:
	rm -f moo
