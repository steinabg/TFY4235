#include <iostream>
#include <string>
#include <omp.h>
// #include <superlu>
#include "stdio.h"
#include "stdlib.h"
#include "list.h"
#include "checkInsideCurveOrg.h"
//#define ARMA_USE_SUPERLU 1;
// #define ARMA_USE_CXX11 1; 
#define ARMA_NO_DEBUG 1;
// #define ARMA_USE_BLAS 1;
#define ARMA_USE_OPENMP 1;
#define ARMA_64BIT_WORD 1;
#include <armadillo>
#include <math.h>
// #include <boost/timer/timer.hpp>
 
 
  
using namespace arma;
using namespace std;


void insertNodesOnSegment(list listOfVertex, int segStart){
// Denne burde flyttes inn i klassen og skrives om til å ikke bruke
	// insert position...

	float x0,y0,x1,y1;
	listOfVertex.get_nodeXY(segStart,&x0,&y0);
	listOfVertex.get_nodeXY(segStart+1,&x1,&y1);
	// printf("x0 = %f, x1 = %f\n", x0,x1);
	float stepx = (x1-x0)/4;
	float stepy = (y1-y0)/4;

	if(stepx){ //True hvis stepx != 0 i.e. horisontal segment
		// Inserting 7 nodes on the line segment
		listOfVertex.insert_position(segStart+1	,x0+stepx   ,y0);
		listOfVertex.insert_position(segStart+2 ,x0+stepx   ,y0+stepx);
		listOfVertex.insert_position(segStart+3 ,x0+2*stepx ,y0+stepx);
		listOfVertex.insert_position(segStart+4 ,x0+2*stepx ,y0);
		listOfVertex.insert_position(segStart+5 ,x0+2*stepx ,y0-stepx);
		listOfVertex.insert_position(segStart+6 ,x0+3*stepx ,y0-stepx);
		listOfVertex.insert_position(segStart+7 ,x0+3*stepx ,y0);
	}
	else{
		// Inserting 7 nodes on the line segment
		listOfVertex.insert_position(segStart+1 ,x0        ,y0+stepy);
		listOfVertex.insert_position(segStart+2 ,x0-stepy  ,y0+stepy);
		listOfVertex.insert_position(segStart+3 ,x0-stepy  ,y0+2*stepy);
		listOfVertex.insert_position(segStart+4 ,x0        ,y0+2*stepy);
		listOfVertex.insert_position(segStart+5 ,x0+stepy  ,y0+2*stepy);
		listOfVertex.insert_position(segStart+6 ,x0+stepy  ,y0+3*stepy);
		listOfVertex.insert_position(segStart+7 ,x0        ,y0+3*stepy);
	}
}

int main(int argc, char **argv){
	int approximationType = 4; // set til 4 eller 9.
	
 
	wall_clock timer;
	timer.tic();
	// printf("Hello, World!\n");
	// printf("The number of arguments is %i and argument 1 in %s\n", argc, argv[1]);
 
  
	int iterations = 4; // Number of generator iterations
	float L = pow(4,iterations); // Length of initial box
	list a; // Initialize linked list to store nodes (x,y) in

	// Create starting geometry createnode(x,y)
	a.createnode(0,0);
	a.createnode(L,0);
	a.createnode(L,-L);
	a.createnode(0,-L);
	a.createnode(0,0); // The curve must be closed for my generator to work
	// boost::timer::auto_cpu_timer t;



	// Iterate over number of orders we want our Koch curve to be
	for (int i = 0; i < iterations; ++i)
	{
		int sides=0;
		sides = a.countSides(); // Find the # of sides in our Koch curve
		for (int i = 0; i < sides; ++i)
		{
			// Apply the generator to all sides
			insertNodesOnSegment(a,i*7+i);
		}

	}
	

	/* insert_spacings(arg) inserts a straight segment of
	 	length arg between corners. */
	a.insert_spacings(2);
	// boost::timer::auto_cpu_timer t1;

	
	// a.writeToFile(); // Coordinates of the list can be written to "coordinates.txt"

	// Define square lattice constant
	/* The lattice constant is given by the distance between two neighbour points*/
	float x0,y0,x1,y1, deltaX, deltaY;
	a.get_nodeXY(0,&x0,&y0);
	a.get_nodeXY(1,&x1,&y1);
	deltaX = x1-x0;
	deltaY = y1-y0;

	float delta = (deltaX > deltaY) ? deltaX : deltaY;
	float minX,minY,maxX,maxY;
	int nodes;
	a.get_boundaries(&minX,&maxX,&minY,&maxY,&nodes);
	printf("minX = %f, minY = %f, maxX = %f, maxY = %f, there are %i nodes in total\n", 
		minX,minY,maxX,maxY,nodes);
	printf("distance between neighbours is %f, which gives %f cols and rows",
		delta,(maxY-minY)/delta+1);

	// Convert the linked list to a matrix
	int rows, columns;
	rows = (maxY-minY)/delta+1;
	columns = (maxX-minX)/delta+1;
	mat GRID(rows,columns, fill::zeros);

	// Fill inn curve nodes in matrix
	//#pragma omp parallel for num_threads(1) // Funker
	for (int i = 0; i < nodes; ++i)
	{
		float x,y;
		a.get_nodeXY(i,&x,&y);
		// We must shift the coordinates, so that they are all positive
		int col = round((x - minX)/delta); // X = 0 -> X = 5
		int row = round((-y + maxY)/delta); 
		GRID(row,col) =  1;
	}

	GRID.save("coordinates.txt", raw_ascii);


	// Print the matrix
	// cout.precision(0);
	// cout.setf(ios::fixed);
	// GRID.raw_print(cout, "GRID:");


	// Test if a point is inside the curve
	// int roww = 7;
	// int coll = 3;
	// cout << "Is (" << roww << " " << coll << ") inside curve? " 
	// << checkInsideCurve(GRID, roww,coll) << endl;
 
   
	// setMatrixInteriorPoints(&GRID);

	// Create sparse matrix to store the number of an interior point in the curve
	sp_mat InteriorNumber(rows,columns);
	int interiorNumberCounter = 0;
	

	// boost::timer::auto_cpu_timer t2;
	// Want to use checkInsideCurve to make interior elements = 2
	int dimension = GRID.n_cols;
	// mat TEMP = GRID; // Hører sammen med OLD SERIAL CODE
 
	
	mat TEMP; 
	TEMP = checkInsideCurveNEW(GRID); // Highly improved performance

   
	GRID = TEMP;
	for (int i = 0; i < dimension; ++i)
	{
		for (int j = 0; j < dimension; ++j)
		{
			if (GRID(i,j) == 2){
				// Gir hvert indre punkt et unikt tall [1,...]
				InteriorNumber(j,i) = ++interiorNumberCounter;
			} 
		}
	} 
 

    /*   OLD SERIAL CODE    */


	// for (int i = 0; i < dimension; ++i)
	// {
	// 	for (int j = 0; j < dimension; ++j)
	// 	{
	// 		if (checkInsideCurve(GRID,i,j)){
	// 			TEMP(i,j) = 2;
	// 			// Gir hvert indre punkt et unikt tall [1,...]
	// 			InteriorNumber(j,i) = ++interiorNumberCounter;

	// 			// Vil ha checkInsideCurve(GRID,TEMP) som sjekker alle 
	// 		} 
	// 	}
	// }
	// GRID = TEMP;


	/*		   NEW PARALLEL CODE   */

	// #pragma omp parallel for num_threads(4) // OBS Må flytte INTERiorno
	// for (int i = 0; i < dimension; ++i)
	// {
	// 	for (int j = 0; j < dimension; ++j)
	// 	{
	// 		if (checkInsideCurve(GRID,i,j)){
	// 			TEMP(i,j) = 2;
	// 			// Gir hvert indre punkt et unikt tall [1,...]
	// 		} 
	// 	}
	// }

	// GRID = TEMP;
	// for (int i = 0; i < dimension; ++i)
	// {
	// 	for (int j = 0; j < dimension; ++j)
	// 	{
	// 		if (GRID(i,j) == 2){
	// 			// Gir hvert indre punkt et unikt tall [1,...]
	// 			InteriorNumber(j,i) = ++interiorNumberCounter;
	// 		} 
	// 	}
	// }

	/*         NEW PARALLEL CODE END    */

	// cout.precision(0);
	// cout.setf(ios::fixed);
	// GRID.raw_print(cout, "GRID:");

	cout << "\nInterior pts = " << interiorNumberCounter << endl;

	// InteriorNumber.print();
	// mat K(InteriorNumber);
	// cout.precision(0);
	// cout.setf(ios::fixed);
	// K.raw_print(cout, "InteriorNumber:");



	// Create coefficient matrix for helmholtz eq
	sp_mat COEFF(interiorNumberCounter,interiorNumberCounter);
	

 
	sp_mat::const_iterator it = InteriorNumber.begin();
	int ii,jj;

	if(approximationType == 4)
	{
		/* FYLL COEFF FOR 5 PKT APPROXIMASJON */
		COEFF.diag() += 4;
		
		for (int i = 0; i < interiorNumberCounter; ++i) // For rader i COEFF
		{

			//Hva er koordinatene (i,j) i GRID, til "-4" i en rad i COEFF?
			ii = it.row();
			jj = it.col();
			// printf("For rad %i, er (ii,jj) = (%i,%i) fra GRID\n",
			// i,ii,jj );
		
			// Sjekk i InteriorNumber om de 4 naboene til interior punktet
			// OBS: Første interiørpunkt skal være 0, derfor har en -1 her.
			int above = InteriorNumber(ii-1,jj)-1;
			int below = InteriorNumber(ii+1,jj)-1;
			int left  = InteriorNumber(ii,jj-1)-1;
			int right = InteriorNumber(ii,jj+1)-1;
			// printf("above = %i\n", above);
			// printf("below = %i\n", below);
			// printf("left = %i\n", left);
			// printf("right = %i\n", right);


			// Set nabopunkter til punkt #i til 1 i coeff
			if (above >= 0) COEFF(i,above) =-1;
			if (below >= 0) COEFF(i,below) =-1;
			if (left  >= 0) COEFF(i,left)  =-1;
			if (right >= 0) COEFF(i,right) =-1;
			
			++it;
		}
	}

	else // 9 pt approximation
	{
		cout << "\nUsing 9 pt approx\n";
		/* TODO:: FYLL COEFF FOR 9 PKT APPROXIMASJON */
		COEFF.diag() += 20;
		// int ii,jj;
		for (int i = 0; i < interiorNumberCounter; ++i) // For rader i COEFF
		{

			//Hva er koordinatene (i,j) i GRID, til "-4" i en rad i COEFF?
			ii = it.row();
			jj = it.col();
			// printf("For rad %i, er (ii,jj) = (%i,%i) fra GRID\n",
			// i,ii,jj );
		
			// Sjekk i InteriorNumber om de 4 naboene til interior punktet
			// OBS: Første interiørpunkt skal være 0, derfor har en -1 her.
			int above = InteriorNumber(ii-1,jj)-1;
			int below = InteriorNumber(ii+1,jj)-1;
			int left  = InteriorNumber(ii,jj-1)-1;
			int right = InteriorNumber(ii,jj+1)-1;
			int rt = InteriorNumber(ii+1,jj-1)-1;
			int lt = InteriorNumber(ii-1,jj-1)-1;
			int rb = InteriorNumber(ii+1,jj+1)-1;
			int lb = InteriorNumber(ii-1,jj+1)-1;
			// printf("above = %i\n", above);
			// printf("below = %i\n", below);
			// printf("left = %i\n", left);
			// printf("right = %i\n", right);


			// Set nabopunkter til punkt #i til 1 i coeff
			if (above >= 0) COEFF(i,above) =-4;
			if (below >= 0) COEFF(i,below) =-4;
			if (left  >= 0) COEFF(i,left)  =-4;
			if (right >= 0) COEFF(i,right) =-4;
			if (rt >= 0) COEFF(i,rt) =-1;
			if (lt >= 0) COEFF(i,lt) =-1;
			if (rb >= 0) COEFF(i,rb) =-1;
			if (lb >= 0) COEFF(i,lb) =-1;
			
			++it;
		}
	}



	cout << "Memory to store COEFF = " << sizeof(COEFF) + COEFF.n_nonzero*sizeof(double) << endl;
	// COEFF.save("COEFF.txt", raw_ascii);
	// COEFF.print();
	// mat Y(COEFF);
	// cout.precision(0);
	// cout.setf(ios::fixed);
	// Y.raw_print(cout, "COEFF:");
	// Y.save("COEFF.txt", raw_ascii);


	double n = timer.toc();
	cout << "\nNumber of seconds: " << n << endl;
	wall_clock timer2;
	timer2.tic();


	cx_vec eigval;
	cx_mat eigvec;

	eigs_gen(eigval,eigvec,COEFF,10,"sm"); //Find eigenvecs and values
	double n2 = timer2.toc();
	cout << "\nNumber of seconds eigs_gen: " << n2 << endl;
	if(approximationType==4) eigval /= (delta*delta);
	else eigval /= 6*(delta*delta);
	// vec eigvalR = real(eigval);
	// eigvalR /= (delta*delta);

	cout << "eigval" << endl;
	eigval.print();
	// cout << "eigvec" << endl;
	// eigvec.print();
	cout << "length = " << eigvec.n_rows << endl;

	


	cube eig1(rows,columns,10, fill::zeros);

	sp_mat::const_iterator it2 = InteriorNumber.begin();
	


	for (int i = 0; i < interiorNumberCounter; ++i) // For rader i COEFF
	{
		ii = it2.row();
		jj = it2.col();
		// printf("i = %d",i);
		for (int j = 0; j < 10; ++j)
		{
			eig1(ii,jj,j) = abs(eigvec(i,j));
			// Skal egentlig være eig1(jj,ii,j), men er rettet i matlab
		}
		// eig1(ii,jj,0) = abs(eigvec(i,3));
		
		++it2; 
	}
	// eig1.print();
	
	/* Print the eigenmodes */
	mat temp;
	std::string name = "eig_";
	for (int i = 0; i < 10; ++i)
	{
		temp = eig1.slice(i);
		temp.save(name + char('0' + i) + ".txt",raw_ascii);
	}
	vec eigtemp = real(eigval);
	eigtemp.save("eigvals.txt",raw_ascii);

	// mat temp;
	// temp = eig1.slice(0);

	// temp.save("eig1.mat", raw_ascii);






	// TODO: Fix so L can be random. The number of matrix elements
	// is incrorrect if I alter L now.








}