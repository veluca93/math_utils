CXXFLAGS:=-O3 -Wall -std=c++17 -flto -g -march=native -I. -fno-exceptions \
	-fopenmp -march=native
LDFLAGS:=-flto -lpthread
CXX:=g++
