CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic -pthread -Icpp/include
PY ?= python3
export PYTHONPATH := python/src

.PHONY: all cpp test test-cpp test-py data reconstruct figures paper clean
all: cpp test

cpp: build/h0audit build/h0audit_tests
build/h0audit: cpp/src/h0audit_cli.cpp $(wildcard cpp/include/h0audit/*.hpp)
	@mkdir -p build && $(CXX) $(CXXFLAGS) $< -o $@
build/h0audit_tests: cpp/tests/test_main.cpp $(wildcard cpp/include/h0audit/*.hpp)
	@mkdir -p build && $(CXX) $(CXXFLAGS) $< -o $@

test: test-cpp test-py
test-cpp: build/h0audit_tests
	./build/h0audit_tests
test-py: build/h0audit
	$(PY) -m unittest discover -s python/tests -v

reconstruct:
	$(PY) scripts/reconstruct_from_paper.py && $(PY) scripts/validate_reconstruction.py

data:            ## needs network; verifies the pinned blob SHA and 30/58 row counts
	$(PY) scripts/fetch_compilation.py
figures:
	$(PY) scripts/make_figures.py
paper: figures
	cd paper && latexmk -pdf -interaction=nonstopmode manuscript.tex
clean:
	rm -rf build && cd paper && latexmk -C >/dev/null 2>&1 || true
