#ifndef BINARY_H
#define BINARY_H
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <utility>   
using namespace std;
typedef pair<int,int> Pair;

//class handeling binary numbers.
//takes an int < 16 and save as a vector of 1 or 0.
//permute value permutes a value-table after a permutation
// of the bits.
class binary{
public:
	binary(){
		bv = vector<unsigned char>(4);
	}
	binary(int n)
	{
		// binary vector, first bit is to the right.
		assert(n < 16 && n >= 0);
		bv = vector<unsigned char>(4);
		bv[3] = (n >> 0) %2;
		bv[2] = (n >> 1) %2;
		bv[1] = (n >> 2) %2;
		bv[0] = (n >> 3) %2;
	}
	void int2bin(int n, int size){
		// binary vector, first bit is to the right.
		assert(n < 16 && n >= 0);
		if(size ==4){
			//cout <<  "111111111111111111111" << endl;
			bv = vector<unsigned char>(4);
			bv[3] = (n >> 0) %2;
			bv[2] = (n >> 1) %2;
			bv[1] = (n >> 2) %2;
			bv[0] = (n >> 3) %2;	
		}

		if(size == 3){
			//cout <<  "22222222222222222222" << endl;
			bv = vector<unsigned char>(3);
			bv[2] = (n >> 0) %2;
			bv[1] = (n >> 1) %2;
			bv[0] = (n >> 2) %2;
		}

		if(size == 2) {
			//cout <<  "33333333333333333333" << endl;
			bv = vector<unsigned char>(2);
			bv[1] = (n >> 0) %2;
			bv[0] = (n >> 1) %2;

		}

	}

	vector<unsigned char> bv;
	void print(){
		cout << (int) bv[3] <<(int) bv[2] << (int) bv[1] << (int) bv[0] << endl;
	}
	// a = 0 , b= 1 is the "normal"
	int permute_value(int a, int b){
		return bv[a] + bv[b]*2;
	}
	int permute_value(int a, int b, int c){
		return bv[a] + bv[b]*2 + bv[c]*4;
	}
	int permute_value(int a, int b, int c, int d){
		return bv[a] + bv[b]*2 + bv[c]*4 + bv[d]*8;
	}

	int permute_value(vector<int> abcd){
		assert(abcd.size() <= 4);
		int value =0;


		if(abcd.size() == 4){
			//cout << "the binary vector: " << (int)bv[0] << (int)bv[1] << (int)bv[2] << (int)bv[3] << endl;
			//cout << "abcd: " << (int)abcd[0] << (int) abcd[1] << (int) abcd[2] << (int) abcd[3] << endl;
			for(int i = 0; i< abcd.size(); i++) value += bv[abcd[i]]* (1<< (abcd.size()-1-i));
			//cout << "value: " << value << endl;
		}

		if(abcd.size() == 3){
			//cout << "the binary vector: " << (int)bv[0] << (int)bv[1] << (int)bv[2] << endl;
			//	cout << "abcd: " << (int)abcd[0] << (int) abcd[1] << (int) abcd[2] << endl;
			for(int i = 0; i< abcd.size(); i++) value += bv[abcd[i]]* (1<< (abcd.size()-1-i));
			//	cout << "value: " << value << endl;
		}
		if(abcd.size() == 2){
			//cout << "the binary vector: " << (int)bv[0] << (int)bv[1] << (int)bv[2] << endl;
			//	cout << "abcd: " << (int)abcd[0] << (int) abcd[1] << (int) abcd[2] << endl;
			for(int i = 0; i< abcd.size(); i++) value += bv[abcd[i]]* (1<< (abcd.size()-1-i));
			//	cout << "value: " << value << endl;
		}

		return value;
	}
};

struct CmpPair
{
	bool operator()(const Pair& a, const Pair& b)
	{ return a.first <  b.first; }
};

void sortingPermutation(vector<int>& values, vector<int>& permutation)
{
	vector<Pair> pairs;
	for (int i = 0; i < (int)values.size(); i++)
		pairs.push_back(Pair(values[i], i));

	std::sort(pairs.begin(), pairs.end(), CmpPair());
	values.clear();
	typedef vector<Pair>::const_iterator I;
	for (I p = pairs.begin(); p != pairs.end(); ++p){
		permutation.push_back(p->second);
		values.push_back(p->first);
	}
}

vector<float> permute_table(int a, int b, const vector<float> & values){
	vector<float> rv(4);
	binary bin;

	for(int i = 0; i< 4 ; i++){
		bin.int2bin(i,2);
		rv.at(i) = values.at(bin.permute_value(a,b));
	}
	return rv;
}
vector<float> permute_table(int a, int b,int c, const vector<float> & values){
	vector<float> rv(9);
	binary bin;

	for(int i = 0; i< 8 ; i++){
		bin.int2bin(i,3);
		rv.at(i) = values.at(bin.permute_value(a,b,c) );
	}
	return rv;
}
vector<float> permute_table(int a, int b,int c,int d, const vector<float> & values){
	vector<float> rv(16);
	binary bin;

	for(int i = 0; i< 16 ; i++){
		bin.int2bin(i,4);
		rv.at(i) = values.at(bin.permute_value(a,b,c,d));
	}
	return rv;
}

void permute_table(const vector<int>& abcd, float* values){
	int n = (1<< abcd.size());
	//cout << "n: "<< n << endl;
	vector<float> rv(n);
	binary bin;


	for(int i = 0; i<n ; i++){
		bin.int2bin(i, abcd.size());   //change here!
		rv.at(bin.permute_value(abcd )) = values[i];

		//	cout << bin.permute_value(abcd ) << ": maps too " << i << endl;
	}
	for(int i = 0; i<n ; i++){
		values[i] = rv[i];

		//	cout << bin.permute_value(abcd ) << ": maps too " << i << endl;
	}
}

void permute_table(const vector<int>& abcd,  vector<float> & values){
	int n = (1<< abcd.size());
	//cout << "n: "<< n << endl;
	vector<float> rv(n);
	binary bin;


	for(int i = 0; i<n ; i++){
		bin.int2bin(i, abcd.size());   //change here!
		rv.at(bin.permute_value(abcd )) = values.at(i);
		//cout << bin.permute_value(abcd ) << ": maps too " << i << endl;
	}
	for(int i = 0; i<n ; i++){
		values[i] = rv[i];

		//	cout << bin.permute_value(abcd ) << ": maps too " << i << endl;
	}
}
#endif