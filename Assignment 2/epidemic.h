#include <iostream>
#include <armadillo>
#include <random>
#include <chrono>
#include "stdlib.h"
#include <string>
#include "stdio.h"




using namespace arma;
using namespace std;


class epidemic
{
private:
	int globalTime;  // Number of timesteps used
	int latticeLength;
	int noMutations; // Maximum number of mutations
	float probOfMutation;
	float ICfractionInfected;
	float p; // Probability of infection
	float q; // Probability of reinfection (q = 0 : Perfect immunization)
	bool stillSickPeople;
	int numberOfInfectedPeopleNow;
	int peopleWhoHasBeenInfected;
	mat saneStatus;	// 0 = frisk, ellers er # = mutasjonsnummer
	sp_mat hasBeenInfected;
	mat carrierStatus;	// 0 = not carrier,  ellers er # = mutasjonsnummer
	cube recoveredStatus;	// Binær kube. (i,j,k) = Individ (i,j) har vært smittet av mutasjon k (1/0)
	mat inbox;	// "Inbox" for individer. (i,j) = mottat mutasjon
	mat timeOfInfection;
	// std::random_device rd;
	mt19937 generator;

public:
	vec returnIndex(int number);
	void sendPathogenToNeighbors(int row, int col);
	void timeStepSpread();
	bool infectiousIndividual(int row, int col);
	bool receivedPathogen(int row, int col);
	void becomeInfected(int row, int col);
	bool hasBeenInfectedBefore(int row, int col, int mutation);
	void setIC();
	void setGenerator(mt19937 gen);
	mt19937 getGenerator();
	int getNumberOfPeopleInfectedNow(); // returns "rho_1"
	int getPeopleWhoHasBeenInfected();  // returns "rho_2"
	void printSaneStatus();
	int runTillExtinct(bool debugPrint); // Runs the simulation until all diseases are extinct. Returns no of timesteps used.
	vec runTimeEvolution(int maxTimeSteps, bool debugPrint);
	epidemic(int latticeLength, float probOfMutation, float ICfractionInfected, float p, float q);
	epidemic(int latticeLength, float probOfMutation, float ICfractionInfected, float p, float q, int noMutations);

	~epidemic();
	
};