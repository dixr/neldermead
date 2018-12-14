all: neldermead

neldermead: main.cpp
	g++ $^ -o $@

clean:
	$(RM) neldermead

rebuild:
	$(MAKE) clean
	$(MAKE)

.PHONY: clean rebuild 
.PRECIOUS: main.cpp utils.hpp
