CXX=g++
CXXFLAGS=-O2 -std=c++11 -Wall -lpthread

all: determinant-thread

determinant: determinant-thread.cpp
								 @echo "compile determinant-thread"
								 $(CXX) $(CXXFLAGS) determinant-thread.cpp -o determinant-thread

clean:
	@echo "clean project"
	rm -rf determinant-thread
	rm -rf *.log
