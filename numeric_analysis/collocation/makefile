SHELL = /bin/bash


SRC = $(wildcard *.cc *.cpp)
EXE = $(basename $(SRC))

default: collocation

collocation: collocation.cc
	g++ collocation.cc -lboost_iostreams -lboost_system -lboost_filesystem -o collocation

