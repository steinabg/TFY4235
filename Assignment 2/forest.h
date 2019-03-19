#include <iostream>
#include <armadillo>
#include <random>
#include <chrono>
#include "stdlib.h"
#include <string>
#include "stdio.h"

//////////////////////////////////////////////////////////////
//	forest er en cube hvor slice 0 inneholder status på alle trær
//	elementene i slice 0 er 0: Ikke tre, 1: Grønt tre, 2: Brennende tre
//	3: dødt tre.
//	slice 1 sier hvilket tidssteg et tre tok fyr.
//////////////////////////////////////////////////////////////


using namespace arma;
using namespace std;


class forest
{
private:
	double probability;
	int forestLength;
	cube FORESTMATRIX;
	float numberoftrees;
	float burnedtrees;
	bool useSlides;
	mt19937 generator;
public:
	void fillForest();	// Populate forest lattice with trees
	void setIC();	// Set the top tow to status = burning.
	bool eastIsGreen(int i, int j);
	bool southIsGreen(int i, int j);
	bool eastIsBurning(int i, int j);
	bool southIsBurning(int i, int j);
	bool myNeighborIsBurning(int i, int j);
	int getEastNeighborIndex(int column);
	void burnForest(int* t);
	void ForestToFile(); // May not work properly here
	void printForestToScreen();
	float getNumberOfTrees();
	float getNumberOfBurnedTrees();
	void setGenerator(mt19937 gen);
	mt19937 getGenerator();
	forest(int L, double p, bool useSlides, mt19937 generator);
	forest(int L, double p, bool useSlides);
	forest();
	~forest();
	
};