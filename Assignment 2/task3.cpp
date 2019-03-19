#include <iostream>
#include <armadillo>
#include <random>
#include <chrono>
#include "stdlib.h"
#include <string>
#include "stdio.h"
#include "epidemic.h"
#include <omp.h>


using namespace arma;
using namespace std;

#define THREADS 4
#define CHUNK 5

// const int latticeLength = 100;	
// const int noMutations = 100;
// const float probOfMutation = 0; // "lambda"
// const float ICfractionInfected = 0.1;
// const float p = 1; // Probability of infection
// const float q = 0; // Probability of reinfection (q = 0 : Perfect immunization)
// std::random_device rd;
// mt19937 generator(rd());
// bool stillSickPeople;
// int numberOfInfectedPeopleNow;
// int peopleWhoHasBeenInfected;


/////////////////////////////////////////////////////////
//		prototypes
void task35();	// COMPLETED
void task36();	// COMPLETED
void task38();	// COMPLETED


int main(int argc, char *argv[]){
	wall_clock timer;
	timer.tic();



	// task35();
	task36();
	// task38();




	double n = timer.toc();
	cout << "\nNumber of seconds: " << n << endl;



}

void task35(){	// COMPLETED
	// Plot rho1 as a function of time for many different p and a few different q. Keep lambda = 0.
	std::random_device rd;
	mt19937 generator1(rd());
	mt19937 generator2(rd());
	mt19937 generator3(rd());
	mt19937 generator4(rd());


	const int latticeLength = 100;	// Needs to be 100 for good plots
	const float probOfMutation = 0; // "lambda"
	const float ICfractionInfected = 0.01;	// Set to 1% as per instructions


	vec pvec = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
	// vec pvec = {1};
	// vec qvec = {0.0,0.5,1};
	vec qvec = {0.0, 0.5, 1.0};

	int maxTimesteps = 1e4;

	int totPeople = latticeLength*latticeLength;


	// Find rho_1's evolution in time
	for (int ii = 0; ii < qvec.n_elem; ++ii)	// For all specified q values
	{
		// mat densityMat(maxTimesteps+1,1, fill::zeros);
		mat densityMat(maxTimesteps+1, pvec.n_elem, fill::zeros);
		// printf("pvec.n_elem = %d\n", pvec.n_elem);
		#pragma omp parallel for schedule(dynamic) num_threads(4)
		for (int jj = 0; jj < pvec.n_elem; ++jj)	// For all specified p values
		{


			// printf("p = %f, q = %f\n", pvec(jj),qvec(ii));
			epidemic e(latticeLength,probOfMutation,ICfractionInfected,pvec(jj),qvec(ii));
			/////////////////////////////////////////////////
			// Give each of the threads a diferent RNG
			if 		(omp_get_thread_num() == 0) e.setGenerator(generator1);
			else if (omp_get_thread_num() == 1) e.setGenerator(generator2);
			else if (omp_get_thread_num() == 2) e.setGenerator(generator3);
			else if (omp_get_thread_num() == 3) e.setGenerator(generator4);
			///////////////////////////////////////////////////
			e.setIC();
 
 
			vec test = e.runTimeEvolution(maxTimesteps,false);
			// test.print("rho_1 = ");
			// densityMat.insert_cols(jj, test );
			densityMat.col(jj) = test;
			// e.printSaneStatus();
			/////////////////////////////////////////////////
			// Reuse the RNG from previous class instance so that the seed is not reset!
			if 		(omp_get_thread_num() == 0) generator1 = e.getGenerator();
			else if (omp_get_thread_num() == 1) generator2 = e.getGenerator();
			else if (omp_get_thread_num() == 2) generator3 = e.getGenerator();
			else if (omp_get_thread_num() == 3) generator4 = e.getGenerator();
			///////////////////////////////////////////////////

		}
		printf("q = %9.2f done!\n", qvec(ii));

		// densityMat.print("densityMat:");

		char arr[sizeof(qvec(ii))];
		sprintf(arr,"%0.2f",qvec(ii));


		char arr2[sizeof(latticeLength)];
		sprintf(arr2,"%d",latticeLength);

		char arr3[sizeof(ICfractionInfected)];
		sprintf(arr3,"%0.2f",ICfractionInfected);

		std::string name = "densityMat_l";
		densityMat.save(name + arr2 + "_q" + arr +"_IC_"+ arr3 +".txt",raw_ascii);



	}

}


void task36(){
	// Find rho1 as a function of time for many different p and one q. Keep lambda = 0.
	//	Set max timestep to a high #. Have pvec "continous" from 0 to 1.
	std::random_device rd;
	mt19937 generator1(rd());
	mt19937 generator2(rd());
	mt19937 generator3(rd());
	mt19937 generator4(rd());

	const int maxTimesteps = 1e4;
	const int latticeLength = 60;	// Needs to be at least 50 for good plots
	const float probOfMutation = 0.2; // "lambda"
	const float ICfractionInfected = 0.01;	// Set to 1% as per instructions
	const int numberOfSimulations = 1; // THIS SHOULD ALWAYS BE 1 HERE

	// vec pvec = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
	vec pvec = linspace(0,1,101);
	vec qvec = {1};
	int saveInterval = ceil(pvec.n_elem*0.1);
	printf("saveInterval = %d\n", saveInterval);
	

	sp_vec checkVector(pvec.n_elem);

	const int noMutations = maxTimesteps + 10;	// We should not get more mutations than there are timesteps
	int totPeople = latticeLength*latticeLength;

	
	cube densityCube(maxTimesteps+1, pvec.n_elem, numberOfSimulations, fill::zeros);

	// Find rho_1's evolution in time
	for (int ii = 0; ii < numberOfSimulations; ++ii)	//
	{
		// mat densityMat(maxTimesteps+1,1, fill::zeros);
		// printf("pvec.n_elem = %d\n", pvec.n_elem);
		#pragma omp parallel for schedule(dynamic) num_threads(4)
		for (int jj = 0; jj < pvec.n_elem; ++jj)	// For all specified p values
		{


			// printf("p = %f, q = %f\n", pvec(jj),qvec(ii));
			epidemic e(latticeLength,probOfMutation,ICfractionInfected,pvec(jj),qvec(0), noMutations);
			/////////////////////////////////////////////////
			// Give each of the threads a diferent RNG
			if 		(omp_get_thread_num() == 0) e.setGenerator(generator1);
			else if (omp_get_thread_num() == 1) e.setGenerator(generator2);
			else if (omp_get_thread_num() == 2) e.setGenerator(generator3);
			else if (omp_get_thread_num() == 3) e.setGenerator(generator4);
			///////////////////////////////////////////////////
			e.setIC();
 
 
			vec test = e.runTimeEvolution(maxTimesteps,false);
			// test.print("rho_1 = ");
			// densityMat.insert_cols(jj, test );
			densityCube.slice(ii).col(jj) = test;
			// e.printSaneStatus();
			/////////////////////////////////////////////////
			// Reuse the RNG from previous class instance so that the seed is not reset!
			if 		(omp_get_thread_num() == 0) generator1 = e.getGenerator();
			else if (omp_get_thread_num() == 1) generator2 = e.getGenerator();
			else if (omp_get_thread_num() == 2) generator3 = e.getGenerator();
			else if (omp_get_thread_num() == 3) generator4 = e.getGenerator();
			///////////////////////////////////////////////////
			// checkVector(ii) = 1;
			// if(omp_get_thread_num() == 0){
			// 	int nNonzero = checkVector.n_nonzero;
			// 	if(!(nNonzero % saveInterval)){
			// 		double percent = (float) nNonzero/pvec.n_elem*4;
			// 		printf("percent done = %0.2f\n", percent);
			// 	}
			
			// }

		}




	}
	mat densityMat = mean(densityCube,2);


	char arr2[sizeof(latticeLength)];
	sprintf(arr2,"%d",latticeLength);

	char arr3[sizeof(qvec(0))];
	sprintf(arr3,"%0.2f",qvec(0));
	char arr4[sizeof(probOfMutation)];
	sprintf(arr4,"%0.2f",probOfMutation);

	std::string name = "rho1Mat_";
	densityMat.save(name + arr2 + "_q"+ arr3 + "lambda" + arr4 +".txt",raw_ascii);


	/////////////////////////////////////////////////////////////
	// const int sampleIntervall = 30;	// # timesteps between each sample
	// const int offset = 600;
	// int noOfSamplePts = (maxTimesteps+1-offset)/sampleIntervall;
	// printf("noOfSamplePts = %d\n", noOfSamplePts);
	// printf("densityMat.n_rows = %d\n", densityMat.n_rows);

	// mat averageMat(noOfSamplePts ,pvec.n_elem, fill::zeros);
	// for (int ii = 0; ii < noOfSamplePts; ++ii)
	// {
	// 	int index = offset + sampleIntervall*ii;
	// 	printf("index = %d\n", index);
	// 	averageMat.row(ii) = densityMat.row(index);
	// }

	// rowvec averageVec = mean(averageMat,0);
	// // printf("averageVec.n_rows = %d\n", averageVec.n_rows);

	// averageVec.save(name + arr2 + "_IC_"+ arr3 +".txt",raw_ascii);



}


void task38(){ // COMPLETED
	// Make q = lambda = 0. Stop the simulation after no more infected individuals.
	// Measure density of individuals who has been infected once, and timesteps used.

	std::random_device rd;
	mt19937 generator1(rd());
	mt19937 generator2(rd());
	mt19937 generator3(rd());
	mt19937 generator4(rd());


	const int latticeLength = 50;			// Tid øker når denne øker
	const int noMutations = 100;
	const float probOfMutation = 0; // "lambda"
	const float ICfractionInfected = 0.1;	// Fin spiss på p = 0.5 når denne er rundt 0.01. Blir med utvasket når den øker.

	// Plot density of infected people as a function of time for different p and a few different q.
	// vec pvec = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
	vec pvec = linspace(0,1,101);
	// vec pvec = {0.3};
	// vec qvec = {0.0,0.5,1};
	vec qvec = {0};
	// int numberOfSimulations = 2000;
	int numberOfSimulations = latticeLength*latticeLength*2;



	int totPeople = latticeLength*latticeLength;
	vec rho_2vec(pvec.n_elem,fill::zeros);
	vec elapsedTimeVec(pvec.n_elem,fill::zeros);
	mat rho_2mat(pvec.n_elem,numberOfSimulations,fill::zeros);
	mat timeMat(pvec.n_elem,numberOfSimulations,fill::zeros);

	omp_lock_t writelock;
	omp_init_lock(&writelock);
	// Find rho_2's evolution in time
	#pragma omp parallel for num_threads(4)
	for (int ii = 0; ii < numberOfSimulations; ++ii)	// For all specified q values
	{


		for (int jj = 0; jj < pvec.n_elem; ++jj)	// For all specified p values
		{


			// printf("p = %f, q = %f\n", pvec(jj),qvec(ii));
			epidemic e(latticeLength,probOfMutation,ICfractionInfected,pvec(jj),0);
			// e.setGenerator(generator);

			/////////////////////////////////////////////////
			// Give each of the threads a diferent RNG
			if 		(omp_get_thread_num() == 0) e.setGenerator(generator1);
			else if (omp_get_thread_num() == 1) e.setGenerator(generator2);
			else if (omp_get_thread_num() == 2) e.setGenerator(generator3);
			else if (omp_get_thread_num() == 3) e.setGenerator(generator4);
			///////////////////////////////////////////////////



			e.setIC();

			// omp_set_lock(&writelock);
			// e.printSaneStatus();
 			

			elapsedTimeVec(jj) = e.runTillExtinct(false);
			rho_2vec(jj) = (float) e.getPeopleWhoHasBeenInfected()/totPeople;
			// omp_unset_lock(&writelock);
			// test.print("rho_1 = ");
			// e.printSaneStatus();
			// generator = e.getGenerator();

			/////////////////////////////////////////////////
			// Reuse the RNG from previous class instance so that the seed is not reset!
			if 		(omp_get_thread_num() == 0) generator1 = e.getGenerator();
			else if (omp_get_thread_num() == 1) generator2 = e.getGenerator();
			else if (omp_get_thread_num() == 2) generator3 = e.getGenerator();
			else if (omp_get_thread_num() == 3) generator4 = e.getGenerator();
			///////////////////////////////////////////////////





		}


		// elapsedTimeVec.print("elapsedTimeVec:");
		// rho_2vec.print("rho_2vec:");

		rho_2mat.col(ii) = rho_2vec;
		// rho_2mat.print("rho_2mat:");

		timeMat.col(ii) = elapsedTimeVec;
		// timeMat.print("timeMat:");
	}
	// rho_2mat.print("rho_2mat:");
	// 	timeMat.print("timeMat:");
	vec averageRho_2 = mean(rho_2mat,1);
	// averageRho_2.print("averageRho_2:");
	vec averageTimes = mean(timeMat,1);
	// averageTimes.print("averageTimes:");


	char arr[sizeof(latticeLength)];
	sprintf(arr,"%d",latticeLength);

	std::string name = "averageRho2_l";
	averageRho_2.save(name + arr + ".txt",raw_ascii);
	std::string name2 = "averageTimes_l";
	averageTimes.save(name2 + arr + ".txt",raw_ascii);


}

