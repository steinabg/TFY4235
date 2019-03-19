#include "forest.h"


forest::~forest(){

}

forest::forest(){
	numberoftrees = 0;
	burnedtrees = 0;
	// seed = chrono::high_resolution_clock::now().time_since_epoch().count();
}


forest::forest(int L, double p, bool useSlides, mt19937 generator): 
forestLength(L), probability(p), useSlides(useSlides), generator(generator)
{
	// forestLength = L;
	// probability = p;
	FORESTMATRIX = cube(forestLength,forestLength,2,fill::zeros);
	numberoftrees = 0;
	burnedtrees = 0;
	// seed = chrono::high_resolution_clock::now().time_since_epoch().count();
}



forest::forest(int L, double p, bool useSlides): 
forestLength(L), probability(p), useSlides(useSlides)
{
	// forestLength = L;
	// probability = p;
	FORESTMATRIX = cube(forestLength,forestLength,2,fill::zeros);
	numberoftrees = 0;
	burnedtrees = 0;
	// seed = chrono::high_resolution_clock::now().time_since_epoch().count();
}


mt19937 forest::getGenerator(){
	return generator;
}


void forest::setGenerator(mt19937 gen){
	generator = gen;
}


float forest::getNumberOfTrees(){
	return numberoftrees;

}
float forest::getNumberOfBurnedTrees(){
	return burnedtrees;
}

void forest::printForestToScreen(){
	FORESTMATRIX.print();
}

void forest::ForestToFile(){ // May not work properly here
	long number = forestLength;
	char buffer [128];
	int ret = snprintf(buffer, sizeof(buffer), "%ld", number);
	char * num_string = buffer; //String terminator is added by snprintf
	// cout << num_string << endl;

	char arr[sizeof(probability)];
	sprintf(arr,"%0.2f",probability);


	std::string name = "forest_size_";
	FORESTMATRIX.slice(0).save(name + num_string + "_p" + arr + ".txt",raw_ascii);



}

void forest::burnForest(int* t){ 
	bool burning = true;
	int time = 1; // So we can find use .n_nonzero to calculate density of burned trees
 
	while(burning)
	{
		burning = false;
		// cout << "time = " << time << endl;

		for (int i = 0; i < forestLength; ++i)
		{
			for (int j = 0; j < forestLength; ++j)
			{
				// printf("i = %d, j = %d\n",i,j);

				if (FORESTMATRIX(i,j,0) == 1) // Green tree
				{
					if (myNeighborIsBurning(i,j)){
						// printf("(%d,%d) WAS set fire to by a NEIGHBOR\n", i,j);
						FORESTMATRIX(i,j,0)++;
						FORESTMATRIX(i,j,1) = time;
						burnedtrees++;
					}

				}

				if (FORESTMATRIX(i,j,0) == 2) //Burning
				{
					// printf("FORESTMATRIX(%d,%d) is burning\n",i,j);
					burning = true;

					if(eastIsGreen(i,j)){
						int jj = getEastNeighborIndex(j);
						FORESTMATRIX(i,jj,0)++;		// Set east neighbor on fire
						FORESTMATRIX(i,jj,1) =time;	// Set the time when the fire started
						burnedtrees++;
						// printf("(%d,%d) sets fire to (%d,%d)\n", i,j,i,jj);
					}
					// cout << "east checked" << endl;
					if(southIsGreen(i,j)){ // Cant be true if i = forestLength-1 so it's ok.
						FORESTMATRIX(i+1,j,0)++;			// Set south neighbor on fire
						FORESTMATRIX(i+1,j,1) =time;		// Set the time when the fire started
						burnedtrees++;
						// printf("(%d,%d) sets fire to (%d,%d)\n", i,j,i+1,j);
					}
					// cout << "south checked" << endl;

					if(FORESTMATRIX(i,j,1) == time-1)
					{
						FORESTMATRIX(i,j,0)++; // If the fire started at the old timestep, put it out.
						// printf("(%d,%d) was put out\n", i,j);
					}
				} 
				// if (FORESTMATRIX(i,j,0) == 1) // Green tree
				// {
				// 	if (myNeighborIsBurning(i,j)){
				// 		// printf("(%d,%d) WAS set fire to by a NEIGHBOR\n", i,j);
				// 		FORESTMATRIX(i,j,0)++;
				// 		FORESTMATRIX(i,j,1) = time;
				// 		burnedtrees++;
				// 	}

				// }

			}
		}
		// system("clear");
		// FORESTMATRIX.print();
		// cin >> a;
		

		time++;
	}
	time -= 2; // Korrigering
	*t = time;
}




void forest::fillForest(){

	if (useSlides)
	{
		mat fromFile;
		fromFile.load("testmatrise.txt", raw_ascii);
		FORESTMATRIX.slice(0) = fromFile;
		uvec q1 = find(FORESTMATRIX.slice(0));
		numberoftrees = q1.n_elem;

	}
	else
	{

		// cube temp = cube1;
		// mt19937 generator(seed);
		// printf("seed = %f\n", seed);
		std::uniform_real_distribution<double> dis(0.0,1.0);

		for (int i = 0; i < forestLength; ++i)
		{
			for (int j = 0; j < forestLength; ++j)
			{
				double value = dis(generator);
				// printf("generator = %f\n", value);
				if (value <= probability)
				{
					FORESTMATRIX(j,i,0) = 1;
					numberoftrees++; // Simultaneuosly count number of trees
				}

			}
			// cout << i;
		}
	}

}

void forest::setIC(){
	// cube temp = cube1;

	double value;
	for (int i = 0; i < forestLength; ++i)
	{
		value = FORESTMATRIX(0,i,0);
		if (value == 1)
		{
			FORESTMATRIX(0,i,0)++;
			FORESTMATRIX(0,i,1)++;
			burnedtrees++;
		}
	}
	// return temp;

}

bool forest::eastIsGreen(int i, int j){
	int jj;
	if(j==forestLength-1) jj = 0;
	else jj = j+1;

	// cout << "jj = " << jj << endl;
	if(FORESTMATRIX(i,jj,0) == 1) return true;
	else return false;
}

bool forest::southIsGreen(int i, int j){
	if( i == forestLength-1 ) return false;
	else if (FORESTMATRIX(i+1,j,0) == 1) return true;
	return false;
}

bool forest::myNeighborIsBurning(int i, int j){
	if(southIsBurning(i,j) || eastIsBurning(i,j)) return true;
	return false;
}

bool forest::southIsBurning(int i, int j){
	if( i == forestLength-1 ) return false;
	else if (FORESTMATRIX(i+1,j,0) == 2) return true;
	return false;
}

bool forest::eastIsBurning(int i, int j){
	int jj;
	if(j==forestLength-1) jj = 0;
	else jj = j+1;

	// cout << "jj = " << jj << endl;
	if(FORESTMATRIX(i,jj,0) == 2) return true;
	else return false;

}

int forest::getEastNeighborIndex(int column){
	if(column==forestLength-1) return 0;
	else return column+1;
}