#include "stdafx.h"
#include "process.h"

void Tau1Tau2TimeSeries(int _Connor,int _HH,int _ML1,int _ML2,int _Inh){
	
	Network neural_network(_Connor,_HH,_ML1,_ML2,_Inh,0);
	neural_network.BuildSquare();

	double g=0.1;
	for(double i=2;i<11;i+=0.2){
		for(double j=1;j<i;j+=0.2){
			neural_network.SetChemical(i,j,g);
			neural_network.OutputTimeSeries();
			printf("i=%lf j=%lf is completed\n",i,j);
		}
	}
			
}