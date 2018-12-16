CXX?=g++

all: neldermead

neldermead: main.cpp
	$(CXX) -pedantic -Wall -Weffc++ -Wextra -Og -g --std=c++11 $^ -o $@

clean:
	$(RM) neldermead

rebuild:
	$(MAKE) clean
	$(MAKE)

.PHONY: clean rebuild
.PRECIOUS: main.cpp utils.hpp
