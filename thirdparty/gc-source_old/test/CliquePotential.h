#pragma once

class CliquePotential
{
public:
	CliquePotential(int clique_size, float* labeling_cost);
	~CliquePotential();

	bool IsSubmodular();

public:
	float* _labelingCosts;
	int _cliqueSize;
	int _numLabeling;
};
