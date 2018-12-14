all: neldermead

neldermead: main.cpp
	g++ -pedantic -Wall -Weffc++ -Wextra -g -std=c++1z $^ -o $@

clean:
	$(RM) neldermead

rebuild:
	$(MAKE) clean
	$(MAKE)

.PHONY: clean rebuild
.PRECIOUS: main.cpp utils.hpp
