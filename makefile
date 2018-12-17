CXX?=g++
CXXFLAGS?=-pedantic -Wall -Weffc++ -Wextra -Og -g --std=c++11

all: neldermead

neldermead: main.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	$(RM) neldermead

rebuild:
	$(MAKE) clean
	$(MAKE)

.PHONY: clean rebuild
.PRECIOUS: main.cpp utils.hpp
