#!/bin/bash

CXX?=g++
CXXFLAGS+= -O3 -std=c++11 -Wall -fopenmp -ffunction-sections -Wl,--gc-sections -flto
DIR=./src/

plsa: $(DIR)plsa.cc $(DIR)docopt.cc
	$(CXX) -o $@ $^ $(CXXFLAGS) 

all: plsa

clean:
	rm -f plsa
