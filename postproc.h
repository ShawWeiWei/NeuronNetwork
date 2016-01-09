#pragma once
#include "Network.h"

typedef struct infoNetwork{
	int Connor_num;
	int HH_num;
	int ML1_num;
	int ML2_num;
    double tau1;
	double tau2;
	double g;
}*pInfoNetwork;
/*Output spatiotemporal (each dot corresponds to the firing of a spike by a neuron)...*/
void OutputSpatioTemporal(const pInfoNetwork pNetwork);
/*Output all spike time...*/
void OutputSpikeTime(const pInfoNetwork pNetwork);
/*Output average frequency of network...*/
void OutputFrequency(const pInfoNetwork pNetwork);
/*Calculate the sum of complex exponents of a array of phases...*/
double CalcSumCompExp(const double*_phase,const int &_S);
/*Printf mean phase coherence...*/
void MeanPhase(const pInfoNetwork pNetwork,const int&_m,const int&_n);
/*Printf synchronous bursting...*/
void SyncBursting(const pInfoNetwork pNetwork);
/*Calculate amplitude of average vector field*/
double AverVecField(const pInfoNetwork pNetwork);
/*Output tau1 tau2 average vector field*/
void Tau1Tau2AverVecField(int _Connor_num,int _HH_num,int _ML1_num, int _ML2_num);
/*Calculate synchronization factor*/
double SyncFactor(const pInfoNetwork pNetwork);
/*Output tau1 tau2 synchronization factor*/
void Tau1Tau2SyncFactor(int _Connor_num,int _HH_num,int _ML1_num, int _ML2_num);