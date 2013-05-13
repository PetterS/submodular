#include "insert_clique2.h"
#include <sstream>

int get_digits2(int number){
	int digits = 0; 
	while (number != 0) {
		number /= 10; digits++; 
	}
	return digits;
}
void add_values(int i, int j, int k, int l, float* val4, float value){
	val4[i] += value;
	val4[j] += value;
	val4[k] += value;
	val4[l] += value;
}
void add_values(int i, int j, float* val4, float value){
	val4[i] += value;
	val4[j] += value;
}

void insert_clique2(int poss, float* val4, float* val2){

	//two variables
	if (get_digits2(poss) == 2){
		switch (poss)
		{
		case 12:
			add_values(0,1,2,3,val4, val2[0]);
			add_values(4,5,6,7,val4, val2[1]);
			add_values(8,9,10,11,val4, val2[2]);
			add_values(12,13,14,15,val4, val2[3]);
			break;

		case 13:
			add_values(0,1,4,5,val4, val2[0]);
			add_values(2,3,6,7,val4, val2[1]);
			add_values(8,9,12,13,val4, val2[2]);
			add_values(10,11,14,15,val4, val2[3]);

			break;
		case 14:
			add_values(0,2,4,6,val4, val2[0]);
			add_values(1,3,5,7,val4, val2[1]);
			add_values(8,10,12,14,val4, val2[2]);
			add_values(9,11,13,15,val4, val2[3]);

			break;
		case 23:
			add_values(0,1,8,9,val4, val2[0]);
			add_values(2,3,10,11,val4, val2[1]);
			add_values(4,5,12,13,val4, val2[2]);
			add_values(6,7,14,15,val4, val2[3]);

			break;
		case 24:
			add_values(0,2,8,10,val4, val2[0]);
			add_values(1,3,9,11,val4, val2[1]);
			add_values(4,6,12,14,val4, val2[2]);
			add_values(5,7,13,15,val4, val2[3]);


			break;
		case 34:
			add_values(0,4,8,12,val4, val2[0]);
			add_values(1,5,9,13,val4, val2[1]);
			add_values(2,6,10,14,val4, val2[2]);
			add_values(3,7,11,15,val4, val2[3]);

			break;
		}
	}

	//three variables (8 values)
	if (get_digits2(poss) == 3){
		switch (poss)
		{
		case 123:
			add_values(0,1,val4, val2[0]);
			add_values(2,3,val4, val2[1]);
			add_values(4,5,val4, val2[2]);
			add_values(6,7,val4, val2[3]);
			add_values(8,9,val4, val2[4]);
			add_values(10,11,val4, val2[5]);
			add_values(12,13,val4, val2[6]);
			add_values(14,15,val4, val2[7]);
			break;

		case 124:
			add_values(0,2,val4, val2[0]);
			add_values(1,3,val4, val2[1]);
			add_values(4,6,val4, val2[2]);
			add_values(5,7,val4, val2[3]);
			add_values(8,10,val4, val2[4]);
			add_values(9,11,val4, val2[5]);
			add_values(12,14,val4, val2[6]);
			add_values(13,15,val4, val2[7]);

			break;
		case 134:
			add_values(0,4,val4, val2[0]);
			add_values(1,5,val4, val2[1]);
			add_values(2,6,val4, val2[2]);
			add_values(3,7,val4, val2[3]);
			add_values(8,12,val4, val2[4]);
			add_values(9,13,val4, val2[5]);
			add_values(10,14,val4, val2[6]);
			add_values(11,15,val4, val2[7]);

			break;
		case 234:
			add_values(0,8,val4, val2[0]);
			add_values(1,9,val4, val2[1]);
			add_values(2,10,val4, val2[2]);
			add_values(3,11,val4, val2[3]);
			add_values(4,12,val4, val2[4]);
			add_values(5,13,val4, val2[5]);
			add_values(6,14,val4, val2[6]);
			add_values(7,15,val4, val2[7]);
			break;
		}
	}

}