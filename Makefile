CXX = g++
CXXFLAGS = -O2 -march=native -Wall -Wextra -std=c++14

SRC = dist_cII.cc            \
      
EXE = $(SRC:.cc=.x)

# eliminate default suffixes
.SUFFIXES:
SUFFIXES =

# just consider our own suffixes
.SUFFIXES: .cc .x

all: $(EXE)

.PHONY: all

%.x: %.cc 
	$(CXX) $< -o $@ $(CXXFLAGS)

format: $(SRC)
	@clang-format -i $^ 2>/dev/null || echo "Please install clang-format to run this command"

.PHONY: format

clean:
	rm -f $(EXE) *~

.PHONY: clean



check: dist_cII.x
	./dist_cII.x ./data/P53dmat/P53dmat_??.bin -o test_output.txt
	@awk '{print $$4}' test_output.txt > tmp_l
	@diff ./reference/labels_aftermerge_ref.txt tmp_l || { echo "ERROR - Labels test failed!"; rm tmp_l; exit 1; }
	@awk '{print $$1,$$2,$$3}' test_output.txt > tmp_dg 
	@diff ./reference/DecisionGraph_ref.txt tmp_dg || { echo "ERROR - Decision graph test failed!"; rm tmp_l; rm tmp_dg; exit 1; }
	@rm -f test_output.txt tmp_l tmp_dg
	@echo "All test passed!"
