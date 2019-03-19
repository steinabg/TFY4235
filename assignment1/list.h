#pragma once
#include <iostream>
#include "stdio.h"
#include "stdlib.h"

struct node
{
	float x;
	float y;
	node *next;
};

class list
{
private:
	node *head, *tail;
public:
	list(){
		head = NULL;
		tail = NULL;
	}
	void writeToFile();
	int countSides();
	void get_boundaries(float* minX, float* maxX, float* minY, float* maxY, int* nodes);
	void get_nodeXY(int pos, float* retX, float* retY);
	void insert_spacings(int spaces);
	// float get_nodeY(int pos);
	void createnode(float x, float y);
	void display();
	void insert_start(float x, float y);
	void insert_position(int pos, float x, float y);
	void delete_first();
	void delete_last();
	void delete_position(int pos);
	
};