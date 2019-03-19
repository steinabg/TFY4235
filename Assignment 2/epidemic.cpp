#include "epidemic.h"

/** epidemic v2.0. Changelog: You can now become reinfected while sick." */

epidemic::~epidemic(){
	
}

epidemic::epidemic(int latticeLength, float probOfMutation, float ICfractionInfected, float p, float q):
latticeLength(latticeLength), probOfMutation(probOfMutation), ICfractionInfected(ICfractionInfected),
p(p), q(q)
{
	// latticeLength = 10;
	noMutations = 100;
	globalTime = 1;
	// probOfMutation = 0;
	// ICfractionInfected = 0.01;
	// p = 1; // Probability of infection
	// q = 0; // Probability of reinfection (q = 0 : Perfect immunization)
	timeOfInfection = mat(latticeLength,latticeLength,fill::zeros);
	hasBeenInfected = sp_mat(latticeLength,latticeLength);
	saneStatus = mat(latticeLength,latticeLength,fill::zeros);	// 0 = frisk, ellers er # = mutasjonsnummer
	carrierStatus = mat(latticeLength,latticeLength,fill::zeros);	// 0 = not carrier,  ellers er # = mutasjonsnummer
	recoveredStatus = cube(latticeLength,latticeLength,noMutations, fill::zeros);	// Binær kube. (i,j,k) = Individ (i,j) har vært smittet av mutasjon k (1/0)
	inbox = mat(latticeLength,latticeLength,fill::zeros);	// "Inbox" for individer. (i,j) = mottat mutasjon
}


epidemic::epidemic(int latticeLength, float probOfMutation, float ICfractionInfected, float p, float q, int noMutations):
latticeLength(latticeLength), probOfMutation(probOfMutation), ICfractionInfected(ICfractionInfected),
p(p), q(q), noMutations(noMutations)
{
	// latticeLength = 10;
	globalTime = 1;
	// noMutations = 100;
	// probOfMutation = 0;
	// ICfractionInfected = 0.01;
	// p = 1; // Probability of infection
	// q = 0; // Probability of reinfection (q = 0 : Perfect immunization)
	timeOfInfection = mat(latticeLength,latticeLength,fill::zeros);
	hasBeenInfected = sp_mat(latticeLength,latticeLength);
	saneStatus = mat(latticeLength,latticeLength,fill::zeros);	// 0 = frisk, ellers er # = mutasjonsnummer
	carrierStatus = mat(latticeLength,latticeLength,fill::zeros);	// 0 = not carrier,  ellers er # = mutasjonsnummer
	recoveredStatus = cube(latticeLength,latticeLength,noMutations, fill::zeros);	// Binær kube. (i,j,k) = Individ (i,j) har vært smittet av mutasjon k (1/0)
	inbox = mat(latticeLength,latticeLength,fill::zeros);	// "Inbox" for individer. (i,j) = mottat mutasjon
}




mt19937 epidemic::getGenerator(){
	return generator;
}

vec epidemic::runTimeEvolution(int maxTimeSteps, bool debugPrint){
	// saneStatus.print("saneStatus before:");
	int timesteps = 0;
	int totPeople = latticeLength*latticeLength;
	vec densityVec(maxTimeSteps+1,fill::zeros);
	densityVec(0) = ICfractionInfected;
	// printf("densityVec.n_elem = %d\n", densityVec.n_elem);

	while(stillSickPeople && timesteps < maxTimeSteps){
		stillSickPeople = false;
		timeStepSpread();
		

		if (debugPrint)
		{
			cout << "\n\n\n";
			cout.precision(0);
			cout.setf(ios::fixed);
			saneStatus.raw_print(cout, "saneStatus:");
			printf("numberOfInfectedPeopleNow = %d\n", numberOfInfectedPeopleNow);
			printf("peopleWhoHasBeenInfected = %d\n", peopleWhoHasBeenInfected);
		}
		timesteps++;
		// printf("densityVec(%d) = %f\n", timesteps,(float) numberOfInfectedPeopleNow/totPeople);
		densityVec(timesteps) = (float) numberOfInfectedPeopleNow/totPeople;


	}
	return densityVec;
}

int epidemic::runTillExtinct(bool debugPrint){
	int timesteps = 0;
	while(stillSickPeople){
		stillSickPeople = false;
		timeStepSpread();

		if (debugPrint)
		{
			cout << "\n\n\n";
			saneStatus.print();
			printf("numberOfInfectedPeopleNow = %d\n", numberOfInfectedPeopleNow);
			printf("peopleWhoHasBeenInfected = %d\n", peopleWhoHasBeenInfected);
		}
		timesteps++;
	}
	return timesteps;
}


void epidemic::setGenerator(mt19937 gen){
	generator = gen;
}

void epidemic::printSaneStatus(){
	saneStatus.print("saneStatus:");
}

int epidemic::getNumberOfPeopleInfectedNow(){
	return numberOfInfectedPeopleNow;
}

int epidemic::getPeopleWhoHasBeenInfected(){
	return peopleWhoHasBeenInfected;
}

void epidemic::setIC(){
	int peopleToInfect = ceil((float) latticeLength*latticeLength*ICfractionInfected);
	std::uniform_int_distribution<> disIC(0, latticeLength*latticeLength-1);
	numberOfInfectedPeopleNow = 0;
	// printf("peopleToInfect = %d\n", peopleToInfect);
	
	while(peopleToInfect != 0)
	{
		/*Infect random people*/
		int rnum = disIC(generator);
		// printf("rnum = %d\n", rnum);
		vec index = returnIndex(rnum);	// Convert from 1D indices to 2D indices
		// index.print("index = ");

		if (!saneStatus(index(0),index(1))) // If not sick
		{
			// printf("saneStatus(%d,%d) = 1\n", index(0), index(1));
			saneStatus(index(0),index(1)) 			= 1; // 1 = mutasjon 0
			carrierStatus(index(0),index(1)) 		= 1;
			recoveredStatus(index(0),index(1),1) 	= 1;
			stillSickPeople = true;
			numberOfInfectedPeopleNow++;
			hasBeenInfected(index(0),index(1)) 	= 1;
			timeOfInfection(index(0),index(1)) = globalTime;


			peopleToInfect--;
		}

	}
	// saneStatus.print("saneStatus:");
	peopleWhoHasBeenInfected = numberOfInfectedPeopleNow;
}

vec epidemic::returnIndex(int number){
	vec index = {0,0};
	index(0) = floor((float) number/latticeLength);
	index(1) = number % (latticeLength);
	return index;
}


void epidemic::timeStepSpread(){
	globalTime++;
	numberOfInfectedPeopleNow = 0;
	// mat saneStatusBefore = saneStatus;

	// Send infection to neighbors
	for (int ii = 0; ii < latticeLength; ++ii)
	{
		for (int jj = 0; jj < latticeLength; ++jj)
		{
			if(infectiousIndividual(ii,jj)){
				sendPathogenToNeighbors(ii,jj);
			}
		}
	}
	// Set status of all who received pathogen to sick (with some probability)
	for (int ii = 0; ii < latticeLength; ++ii)
	{
		for (int jj = 0; jj < latticeLength; ++jj)
		{
			if(receivedPathogen(ii,jj)){ // Is a person sane and has received some pathogen?
				becomeInfected(ii,jj);
				// inbox(ii,jj) = 0;
				// printf("(%d,%d) has received a pathogen\n",ii,jj );
			}
			if(timeOfInfection(ii,jj) == globalTime-1 && timeOfInfection(ii,jj) > 0){
				saneStatus(ii,jj) = 0;
				carrierStatus(ii,jj) = 0;
				// numberOfInfectedPeopleNow--;
			}
			if (saneStatus(ii,jj) != 0)
			{
				numberOfInfectedPeopleNow++;
			}
			
		}
	}
	peopleWhoHasBeenInfected = hasBeenInfected.n_nonzero;

	// Reset the inbox
	inbox = zeros(latticeLength,latticeLength);

}

void epidemic::becomeInfected(int row, int col){
	std::uniform_real_distribution<double> dis(0,1); // For probability of mutation

	int receivedPathogen = inbox(row,col);

	// We want the pathogen to have a chance of mutating before infecting the person
	int pathogenToBeInfectedWith = (dis(generator) > probOfMutation) ? receivedPathogen : (receivedPathogen+1);
	// int pathogenToBeInfectedWith;
	// float rng = dis(generator);
	// if (rng <= probOfMutation)
	// {
	// 	pathogenToBeInfectedWith = receivedPathogen+1;
	// 	printf("mutated!\n");
	// }
	// else pathogenToBeInfectedWith = receivedPathogen;
	// printf("rng = %f and probOfMutation = %f\n", rng, probOfMutation );

	if(pathogenToBeInfectedWith > recoveredStatus.n_slices){ 
		pathogenToBeInfectedWith = 1; // Reset mutation if it exceeds the maximum no of mutations
	}
	// printf("pathogenToBeInfectedWith = %d\n", pathogenToBeInfectedWith);


	if(hasBeenInfectedBefore(row,col,receivedPathogen)) // Check if the person has been infected before
	{
		if (dis(generator) <= p*q) // Become infected with some probability p*q
		{
			saneStatus(row,col) = pathogenToBeInfectedWith;
			carrierStatus(row,col) = pathogenToBeInfectedWith;
			recoveredStatus(row,col,pathogenToBeInfectedWith) = 1;
			// numberOfInfectedPeopleNow++;
			stillSickPeople = true;
			hasBeenInfected(row,col) 	= 1;
			timeOfInfection(row,col) = globalTime;
		}
	}
	else
		if (dis(generator) <= p)	// Become infected with some probability p
		{
			saneStatus(row,col) = pathogenToBeInfectedWith;
			carrierStatus(row,col) = pathogenToBeInfectedWith;
			recoveredStatus(row,col,pathogenToBeInfectedWith) = 1;
			// numberOfInfectedPeopleNow++;
			stillSickPeople = true;
			hasBeenInfected(row,col) 	= 1;
			timeOfInfection(row,col) = globalTime;
		}

}

bool epidemic::hasBeenInfectedBefore(int row, int col, int mutation){
	if(recoveredStatus(row,col,mutation)) return true;
	return false;
}


bool epidemic::receivedPathogen(int row, int col){
	if(inbox(row,col) != 0) return true;
	return false;
}

bool epidemic::infectiousIndividual(int row, int col){
	if(saneStatus(row,col) != 0 || carrierStatus(row,col) != 0) return true;
	return false;
}

void epidemic::sendPathogenToNeighbors(int row, int col){ //Sends the mutation number to the "inbox" of the neighbors
	int myMutation = (saneStatus(row,col) !=0) ? saneStatus(row,col) : carrierStatus(row,col);

	///////////////////////////////////////////////
	// Get correct indices according to PBC //
	int above = (row > 0) ? row-1 : latticeLength-1;
	int below = (row < latticeLength-1) ? row+1 : 0;
	int east  = (col < latticeLength-1) ? col+1 : 0;
	int west  = (col > 0) ? col-1 : latticeLength-1;
	////////////////////////////////////////////////

	////////////////////////////////////////////////
	// Assume neighbors are the 4 nearest squares
	// inbox(above,west) = myMutation;
	inbox(above,col) = myMutation;
	// inbox(above,east) = myMutation;
	inbox(row,west) = myMutation;
	inbox(row,east) = myMutation;
	// inbox(below,west) = myMutation;
	inbox(below,col) = myMutation;
	// inbox(below,east) = myMutation;
}

