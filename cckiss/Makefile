cckiss/cckiss : cckiss/cckiss.cc
	$(CXX) -std=c++14 -O2 cckiss/cckiss.cc -o cckiss/cckiss -Wall -Wextra -Wno-constant-logical-operand -Wno-unused-function

.cckiss.PHONY :
	@echo -n ""

cckiss/%.c.s : .cckiss.PHONY cckiss/cckiss
	cckiss/cckiss $@ $(CC) .cckiss.CPPFLAGS $(CPPFLAGS) .cckiss.CXXFLAGS $(CFLAGS)

cckiss/%.c.o : .cckiss.PHONY cckiss/cckiss
	cckiss/cckiss $@ $(CC) .cckiss.CPPFLAGS $(CPPFLAGS) .cckiss.CXXFLAGS $(CFLAGS)

cckiss/%.s : .cckiss.PHONY cckiss/cckiss
	cckiss/cckiss $@ $(CXX) .cckiss.CPPFLAGS $(CPPFLAGS) .cckiss.CXXFLAGS $(CXXFLAGS)

cckiss/%.o : .cckiss.PHONY cckiss/cckiss
	cckiss/cckiss $@ $(CXX) .cckiss.CPPFLAGS $(CPPFLAGS) .cckiss.CXXFLAGS $(CXXFLAGS)

cckiss/Makefile:
	@echo -n ""
