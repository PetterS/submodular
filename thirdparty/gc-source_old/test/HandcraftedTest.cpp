#include "StdAfx.h"
#include "HandcraftedTest.h"
#include "prgc/PRGC.h"
#include "apgc/APGC.h"
#include "test/CliquePotential.h"
using namespace std;

void HandcraftedTest::Test()
{
	//GadgetGraph* ncgraph = new PRGC(4,1,4,1);
	GadgetGraph* ncgraph = new APGC(4,1,4,1);

	int vars[4];
	float labeling_cost[16];

	labeling_cost[0] = 70;
	labeling_cost[1] = 141;
	labeling_cost[2] = 141;
	labeling_cost[3] = 71;
	labeling_cost[4] = 141;
	labeling_cost[5] = 141;
	labeling_cost[6] = 200;
	labeling_cost[7] = 141;
	labeling_cost[8] = 141;
	labeling_cost[9] = 200;
	labeling_cost[10] = 141;
	labeling_cost[11] = 141;
	labeling_cost[12] = 141;
	labeling_cost[13] = 141;
	labeling_cost[14] = 141;
	labeling_cost[15] = 40;

	CliquePotential* potential = new CliquePotential(4, labeling_cost);
	printf("\nchecking if potential is sub-modular...\n");
	bool submodular = potential->IsSubmodular();
	if(submodular){
		PRINTF("potential submodular \n\n");
	}else{
		PRINTF("potential NOT submodular \n\n");
	}
	delete potential;

	vars[0] = 0; vars[1] = 1; vars[2] = 2; vars[3] = 3;
	ncgraph->AddHigherTerm(vars, labeling_cost);

	ncgraph->AddUnaryTerm(0, 0, 100);
	ncgraph->AddUnaryTerm(1, 0, 100);
	ncgraph->AddUnaryTerm(2, 200, 0);
	ncgraph->AddUnaryTerm(3, 0, 50);

	double flow_value = ncgraph->FindMaxFlow();

	printf("\nfinal labels \t");
	for(int i=0; i<4; ++i)
		printf("%d ", ncgraph->GetLabel(i));
	printf("\n");
}
