#include "list.h"
#include <iostream>
#include <fstream>
#include "stdio.h"
#include "stdlib.h"

using namespace std;

void list::insert_spacings(int spaces){


	node* tempStart= new node;
	node* tempFin= new node;
	tempStart = head;
	tempFin = head->next;
	
	// float x0 = tempStart->x;
	// float y0 = tempStart->x;
	// float x1 = tempFin->x;
	// float y1 = tempFin->y;

	// printf("x0 = %f, x1 = %f\n", x0,x1);




	while(tempFin!=NULL){
		
		// Find the coordinates of the nearby corner-nodes

		float x0 = tempStart->x;
		float y0 = tempStart->y;
		float x1 = tempFin->x;
		float y1 = tempFin->y;

		float deltax = (x1-x0);
		float deltay = (y1-y0);
		float delta = ((deltax == 0) ? deltay : deltax);


		float interMedx,interMedy; //Storage for coordinates of new pts
		// printf("tempStart is (%f,%f) and tempFin is (%f,%f)\n", x0,y0,x1,y1);
	
			
		// Maybe assign the difference to a variable to save computation
		if(x1-x0 != 0) // I.e. horisontal segment
		{
			interMedx = delta/(spaces+1);
			interMedy = 0;

		}
		else // Vertical segment
		{
			interMedx = 0;
			interMedy = delta/(spaces+1);
		}
		for (int i = 0; i < spaces; ++i) // For no. of pts to be inserted
		{

			node* pre=new node;
			node* cur=new node;
			node* temp= new node;
			cur = tempStart;
			for (int j = 0; j < i+1; ++j)
			{
				pre = cur;
				cur=cur->next;
				// printf("switching cur j = %i, i = %i\n", j,i);
			}
			temp->x = x0 + (i+1) * interMedx;
			temp->y = y0 + (i+1) * interMedy;
			// printf("Inserting point at (%f,%f), between (%f,%f) and (%f,%f)\n", 
				// x0 + (i+1) * interMedx,y0 + (i+1) * interMedy, x0,y0,x1,y1);
			pre->next = temp;
			temp->next = cur;
			// display();
		}

		// for (int i = 0; i < spaces+1; ++i) 
		// {
		// 	//Skip the inserted nodes
		// 	tempStart = tempStart->next;
		// 	// printf("tempstart is switching\n");
		// }
		tempStart = tempFin;
		
		tempFin = tempFin->next;
	}



}


void list::get_boundaries(float* minX, float* maxX, float* minY, float* maxY, int* nodes){
	node* temp= new node;
	temp = head;
	*nodes = 0;
	while(temp!=NULL){
		if(*minX > temp->x) *minX = temp->x;
		if(*minY > temp->y) *minY = temp->y;
		if(*maxX < temp->x) *maxX = temp->x;
		if(*maxY < temp->y) *maxY = temp->y;
		temp = temp->next;
		*nodes = *nodes +1;
	}

}


void list::writeToFile(){
	ofstream myfile;
	myfile.open("coordinates.txt");

	node* temp= new node;
	temp = head;
	while(temp!=NULL){
		myfile << temp->x << "\t" << temp->y << "\n";
		temp = temp->next;
	}
	myfile.close();

}

int list::countSides(){
	node* temp= new node;
	temp = head;
	int counter=0;
	while(temp!=NULL){
		temp = temp->next;
		counter++;
	}
	return --counter;
}

void list::createnode(float x, float y){
	node *temp = new node;
	temp->x=x;
	temp->y=y;
	temp->next=NULL;
	if(head==NULL){
		head=temp;
		tail=temp;
		temp=NULL;
	}
	else{
		tail->next = temp;
		tail = temp;
	}



}

void list::display(){
	node* temp= new node;
	temp = head;
	while(temp!=NULL){
		cout << "x = " << temp->x << ", y = " << temp->y << "\t";
		temp = temp->next;
	}
	cout << endl;
}

void list::insert_start(float x, float y){
	node* temp = new node;
	temp->x = x;
	temp->y = y;
	temp->next = head;
	head = temp;

}
void list::get_nodeXY(int pos,float* retX, float* retY){
	node* cur=new node;
	cur = head;
	for (int i = 0; i < pos; ++i)
	{
		cur=cur->next;
		
	}
	*retX = cur->x;
	*retY = cur->y;


}
// float list::get_nodeY(int pos){
// 	node* cur=new node;
// 	cur = head;
// 	for (int i = 0; i < pos; ++i)
// 	{
// 		cur=cur->next;
// 	}
// 	return cur->x;
// }

void list::insert_position(int pos, float x, float y){
	node* pre=new node;
	node* cur=new node;
	node* temp= new node;
	cur = head;
	for (int i = 0; i < pos; ++i)
	{
		pre = cur;
		cur=cur->next;
	}
	temp->x = x;
	temp->y = y;
	pre->next = temp;
	temp->next = cur;

}

void list::delete_first(){
	node* temp = new node;
	temp = head;
	head = head->next;
	delete temp;

}

void list::delete_last(){
	node* current = new node;
	node* previous = new node;
	current = head;
	while(current->next!=NULL){
		previous=current;
		current=current->next;
	}
	tail=previous;
	previous->next=NULL;
	delete current;

}

void list::delete_position(int pos){
	node* current = new node;
	node* previous = new node;
	current = head;
	for (int i = 0; i < pos; ++i)
	{
		previous=current;
		current=current->next;
	}
	previous->next=current->next;
}