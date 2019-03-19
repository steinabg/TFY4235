#include <iostream>
#include <armadillo>
#include <random>
#include <chrono>
#include "stdlib.h"
#include <string>
#include "stdio.h"
#include "forest.h"
#include <omp.h>


using namespace arma;
using namespace std;

#define THREADS 4
#define CHUNK 5



void task21();
void task23();
void task24();
void task25();


int main(int argc, char *argv[]){
	wall_clock timer;
	timer.tic();



	// task21();
	task23();
	// task24();
	// task25();


		//////////////////////////////////////////////////////////
		// OBJECT ORIENTED
		// prototype: forest(int L, double p, bool useSlidesExample, mt19937 generator);
		// forest test(forestLength,probability,0, generator1);
		// test.fillForest();
		// test.setIC();
		// int time;
		// test.burnForest(&time);
		// test.printForestToScreen();
		// float totalTrees = test.getNumberOfTrees();
		// float burnedTrees = test.getNumberOfBurnedTrees();
		// test.ForestToFile();

		// float density = burnedTrees/totalTrees;
		// cout << "Timesteps = " << time << endl;
		// printf("Density = %f\n", density);
		// cout << "burnedTrees = " << burnedTrees;
		// cout << "\ntotaltrees = " << totalTrees << endl;


  



	double n = timer.toc();
	cout << "\nNumber of seconds: " << n << endl;



}

void task21(){
	int forestLength = 100;
	float p = 0.8;
	long seed = chrono::high_resolution_clock::now().time_since_epoch().count();
	mt19937 generator1(seed); // Use different generators for each thread each with a different seed 
	forest f(forestLength,p,0,generator1);
	f.fillForest();
	f.ForestToFile();

}

void task23(){
	int forestLength = 10;
	vec pvec = {0.2,0.4,0.6,0.8};
	for (int ii = 0; ii < pvec.n_elem; ++ii)
	{
		long seed = chrono::high_resolution_clock::now().time_since_epoch().count();
		mt19937 generator1(seed); // Use different generators for each thread each with a different seed 
		forest f(forestLength,pvec(ii),0,generator1);
		f.fillForest();
		f.setIC();
		int time;
		f.burnForest(&time);
		f.printForestToScreen();
		cout << "\n\n\n";
		f.ForestToFile();
	}
}

void task24(){
	std::random_device rd;
	
	vec forestLengthVec = {70};

	for (int kk = 0; kk < forestLengthVec.n_elem; ++kk) // Loop for performing same operation for several forest lengths
	{
		// Declare variables
		double probability,totalTrees, burnedTrees; 
		int ii,jj, time;


		int forestLength = forestLengthVec(kk);
		long seed = chrono::high_resolution_clock::now().time_since_epoch().count();
		mt19937 generator1(seed+rd()); // Use different generators for each thread each with a different seed 
		mt19937 generator2(seed+rd()); // Use different generators for each thread each with a different seed
		mt19937 generator3(seed+rd()); // Use different generators for each thread each with a different seed
		mt19937 generator4(seed+rd()); // Use different generators for each thread each with a different seed


		
		int MonteCarloLength = forestLength*forestLength*2; // Regular no. of steps in MC integration
		// int MonteCarloLength = 100;

		int saveInterval = ceil((float) MonteCarloLength*0.1);
		vec pspace = linspace<vec>(0,1,101); // All probabilities
		// vec pspace = linspace<vec>(0.55,0.65,11);	// Probability vector
		int lengthXspace = pspace.n_elem;	// Number of probabilities

		// printf("lengthXspace = %d\n", lengthXspace);

		mat densities(lengthXspace,MonteCarloLength,fill::zeros);
		mat timesteps(lengthXspace, MonteCarloLength,fill::zeros);


		for (ii = 0; ii < MonteCarloLength; ++ii) // Antall ganger sannsynligheten skal varieres
		{
			#pragma omp parallel for schedule(dynamic, CHUNK) num_threads(THREADS) \
		private(time, probability, totalTrees, burnedTrees) shared(densities,timesteps, forestLength)
			for (jj = 0; jj < lengthXspace; ++jj) // For alle sannsynligheter
			{
				// probability = (double) jj/100;
				probability = pspace(jj);
				// if(probability ==101) cout << "yes\n";
				// printf("%f lagres på index jj = %d\n", probability,jj);
				// double probability = 0.9; // for testing av thread safe random generator
				forest test(forestLength,probability,0);

				/////////////////////////////////////////////////
				// Give each of the threads a diferent RNG
				if 		(omp_get_thread_num() == 0) test.setGenerator(generator1);
				else if (omp_get_thread_num() == 1) test.setGenerator(generator2);
				else if (omp_get_thread_num() == 2) test.setGenerator(generator3);
				else if (omp_get_thread_num() == 3) test.setGenerator(generator4);
				///////////////////////////////////////////////////

				test.fillForest();
				test.setIC(); 
				// int time;
				test.burnForest(&time);
				totalTrees = test.getNumberOfTrees();
				burnedTrees = test.getNumberOfBurnedTrees();
				densities(jj,ii) = burnedTrees/totalTrees; // OBS: false sharing??? SJEKK
				timesteps(jj,ii) = time;

				/////////////////////////////////////////////////
				// Reuse the RNG from previous class instance so that the seed is not reset!
				if 		(omp_get_thread_num() == 0) generator1 = test.getGenerator();
				else if (omp_get_thread_num() == 1) generator2 = test.getGenerator();
				else if (omp_get_thread_num() == 2) generator3 = test.getGenerator();
				else if (omp_get_thread_num() == 3) generator4 = test.getGenerator();
				///////////////////////////////////////////////////
			}

			if(!(ii % saveInterval)){ // Backup the work done

				double percent = (double) ii/(MonteCarloLength);
				char arr[sizeof(probability)];
				sprintf(arr,"%0.2f",percent);
				printf("Percentage done: %0.2f\n",percent);
				std::string name = "densities_completed_";
				std::string name2 = "timesteps_completed_";
				densities.save(name + arr + ".txt",raw_ascii);
				timesteps.save(name2 + arr + ".txt",raw_ascii);
				
			}

		}
		// Save the work to disk
		vec aveDensity = mean(densities,1);
		vec aveTime = mean(timesteps,1);

		char arr[sizeof(forestLength)];
		sprintf(arr,"%d",forestLength);

		std::string name =  "AverageDensities_length_";
		std::string name2 = "AverageTimesteps_length_";
		aveDensity.save(name +  arr + ".txt",raw_ascii);
		aveTime.save(name2 + arr + ".txt",raw_ascii);

		cout << "\n\n";
	}


}


void task25(){
	std::random_device rd;
	
	vec forestLengthVec = {20,30,40,50,60,70};

	for (int kk = 0; kk < forestLengthVec.n_elem; ++kk)
	{
		double probability,totalTrees, burnedTrees; 
		int ii,jj, time;


		int forestLength = forestLengthVec(kk);
		long seed = chrono::high_resolution_clock::now().time_since_epoch().count();
		mt19937 generator1(seed+rd()); // Use different generators for each thread each with a different seed 
		mt19937 generator2(seed+rd()); // Use different generators for each thread each with a different seed
		mt19937 generator3(seed+rd()); // Use different generators for each thread each with a different seed
		mt19937 generator4(seed+rd()); // Use different generators for each thread each with a different seed


		
		int MonteCarloLength = forestLength*forestLength*2; // Regular no. of steps in MC integration
		// int MonteCarloLength = 100;

		int saveInterval = ceil((float) MonteCarloLength*0.1);
		// vec pspace = linspace<vec>(0,1,101); // All probabilities
		vec pspace = linspace<vec>(0.55,0.65,11);	// Probability vector
		int lengthXspace = pspace.n_elem;	// Number of probabilities

		// printf("lengthXspace = %d\n", lengthXspace);

		mat densities(lengthXspace,MonteCarloLength,fill::zeros);
		mat timesteps(lengthXspace, MonteCarloLength,fill::zeros);


		for (ii = 0; ii < MonteCarloLength; ++ii) // Antall ganger sannsynligheten skal varieres
		{
			#pragma omp parallel for schedule(dynamic, CHUNK) num_threads(THREADS) \
		private(time, probability, totalTrees, burnedTrees) shared(densities,timesteps, forestLength)
			for (jj = 0; jj < lengthXspace; ++jj) // For alle sannsynligheter
			{
				// probability = (double) jj/100;
				probability = pspace(jj);
				// if(probability ==101) cout << "yes\n";
				// printf("%f lagres på index jj = %d\n", probability,jj);
				// double probability = 0.9; // for testing av thread safe random generator
				forest test(forestLength,probability,0);

				/////////////////////////////////////////////////
				// Give each of the threads a diferent RNG
				if 		(omp_get_thread_num() == 0) test.setGenerator(generator1);
				else if (omp_get_thread_num() == 1) test.setGenerator(generator2);
				else if (omp_get_thread_num() == 2) test.setGenerator(generator3);
				else if (omp_get_thread_num() == 3) test.setGenerator(generator4);
				///////////////////////////////////////////////////

				test.fillForest();
				test.setIC(); 
				// int time;
				test.burnForest(&time);
				totalTrees = test.getNumberOfTrees();
				burnedTrees = test.getNumberOfBurnedTrees();
				densities(jj,ii) = burnedTrees/totalTrees; // OBS: false sharing??? SJEKK
				timesteps(jj,ii) = time;

				/////////////////////////////////////////////////
				// Reuse the RNG from previous class instance so that the seed is not reset!
				if 		(omp_get_thread_num() == 0) generator1 = test.getGenerator();
				else if (omp_get_thread_num() == 1) generator2 = test.getGenerator();
				else if (omp_get_thread_num() == 2) generator3 = test.getGenerator();
				else if (omp_get_thread_num() == 3) generator4 = test.getGenerator();
				///////////////////////////////////////////////////
			}

			if(!(ii % saveInterval)){

				double percent = (double) ii/(MonteCarloLength);
				char arr[sizeof(probability)];
				sprintf(arr,"%0.2f",percent);
				printf("Percentage done: %0.2f\n",percent);
				std::string name = "densities_completed_";
				std::string name2 = "timesteps_completed_";
				densities.save(name + arr + ".txt",raw_ascii);
				timesteps.save(name2 + arr + ".txt",raw_ascii);
				
			}

		}

		vec aveDensity = mean(densities,1);
		vec aveTime = mean(timesteps,1);

		char arr[sizeof(forestLength)];
		sprintf(arr,"%d",forestLength);

		std::string name =  "AverageDensities_length_";
		std::string name2 = "AverageTimesteps_length_";
		aveDensity.save(name +  arr + ".txt",raw_ascii);
		aveTime.save(name2 + arr + ".txt",raw_ascii);

		cout << "\n\n";
	}


}