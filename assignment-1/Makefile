CXX=g++
CXXFLAGS=-O2 -std=c++11 -Wall -lpthread

all: determinant_thread

determinant_thread: determinant_thread.cpp
								 @echo "compile determinant_thread"
								 $(CXX) $(CXXFLAGS) determinant_thread.cpp -o determinant_thread

clean:
	@echo "clean project"
	rm -rf determinant_thread
	rm -rf *.log
