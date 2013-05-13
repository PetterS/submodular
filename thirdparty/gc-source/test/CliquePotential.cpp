#include "stdafx.h"
#include "CliquePotential.h"
	
CliquePotential::CliquePotential(int clique_size, float* labeling_cost)
:_cliqueSize(clique_size)
{
	_numLabeling = cvRound(pow(2.0, clique_size));
	_labelingCosts = new float[_numLabeling];
	for(int i=0; i<_numLabeling; ++i)
		_labelingCosts[i] = labeling_cost[i];
}

CliquePotential::~CliquePotential(){
	delete[] _labelingCosts;
}

bool CliquePotential::IsSubmodular()
{
	MY_ASSERT(_cliqueSize <= 16);

	float e_x, e_y, e_or, e_and;
	int num_subsets = cvRound(pow(2.0f,_cliqueSize));
	for(int i=0; i<num_subsets; ++i)
	{
		e_x = _labelingCosts[i];
		for(int j=i+1; j<num_subsets; ++j)
		{
			e_y = _labelingCosts[j];
			e_or = _labelingCosts[i | j];
			e_and = _labelingCosts[i & j];
			if((e_x + e_y) < (e_or + e_and)){
				PRINTF("x = %s\n", D_BINARY(i));
				PRINTF("y = %s\n", D_BINARY(j));
				//std::cout << "x = " << i << std::endl;
				//std::cout << "y = " << j << std::endl;
				std::cout << "e_x = " << e_x << ", e_y = " << e_y << ", e_or = " << e_or << ", e_and = " << e_and << std::endl;
				return false;
			}
		}
	}
	return true;
}
