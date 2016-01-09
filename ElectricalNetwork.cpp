#include "stdafx.h"
#include "ElectricalNetwork.h"

ElectricalSNetwork::ElectricalSNetwork():SparseNetwork(){
}

ElectricalSNetwork::ElectricalSNetwork(int _Connor,int _HH,int _ML1,int _ML2):SparseNetwork(_Connor,_HH,_ML1,_ML2,0,0){
}



void ElectricalSNetwork::SetElectrical(const double _gc_e){
	m_coupleintensity=_gc_e;
}

void ElectricalSNetwork::UpdateCouple(){

	memset(sm_gCouple,0,sizeof(double)*m_nNeuron);
	int Total_Number,index;
	double fCouple,fV;
	for(int i=0;i<m_nNeuron;++i){
		Total_Number=m_gCoupleNum[i];
		fCouple=0.0;
		fV=m_gNeuron[i]->V;
		for(int j=0;j<Total_Number;++j){
			index=m_gCoupleIdx[i][j];
			fCouple+=(m_gNeuron[index]->V-fV);
		}
		sm_gCouple[i]=m_coupleintensity*fCouple;
	}
}

void ElectricalSNetwork::contextualDirectory(char *direct){
//	char dir[40];
	sprintf_s(direct,100,"F:\\output\\Electrical\\%s\\Connor%dHH%dML1_%dML2_%d",m_sName,m_nConnor,m_nHH,m_nML1,m_nML2);
	_mkdir("F:\\output\\Electrical");
	char direct1[50];
	sprintf_s(direct1,"F:\\output\\Electrical\\%s",m_sName);
	_mkdir(direct1);
	_mkdir(direct);
}

void ElectricalSNetwork::SpecifyFilename(char *filename){
	sprintf_s(filename,40,"gc=%.5lf",m_coupleintensity);
}