#include "stdafx.h"
#include "Network.h"
#include "MorrisLecarModel.h"
#include "HHModel.h"
#include "ConnorModel.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
using namespace std;
//静态成员变量初始化
double Network::sm_dt=0.1;
double *Network::sm_gCouple=NULL;
double *Network::sm_gNoise=NULL;


Network::Network():m_nConnor(0),m_nHH(0),m_nML1(100),m_nML2(0),m_nML2_Inh(0),m_nML2_Inh_Isolated(0),m_nNeuron(100){

	
	/*dynamically allocate memory for coupling matrix*/
	m_gProperty=new u8[m_nNeuron];
	memset(m_gProperty,0,sizeof(u8)*m_nNeuron);
	m_sName=new char[40];
	memset(m_sName,0,sizeof(char)*40);
	if(sm_gCouple==NULL){
		sm_gCouple=new double[m_nNeuron];
		memset(sm_gCouple,0,sizeof(double)*m_nNeuron);
	}
	if(sm_gNoise==NULL){
		sm_gNoise=new double[m_nNeuron];
		memset(sm_gNoise,0,sizeof(double)*m_nNeuron);
	}
//	printf("CoupleVector is %p\n",CoupleVector);
	/*dynamically allocate memory for instances of neuron model*/
	m_gNeuron=new neuron*[m_nNeuron];
	for(int i=0;i<m_nNeuron;++i){
		m_gNeuron[i]=new MorrisLecar(i,1,1);
		m_gProperty[i]=NF_CLASS_ONE|NF_RESTING|NF_ISCOUPLED|NF_EXCITATION;
	}
}

Network::Network(const int _Connor,const int _HH,const int _ML1,const int _ML2,const int _INHIB,const int _INHIB_NOCOUPLING):m_nConnor(_Connor),m_nHH(_HH)
	,m_nML1(_ML1),m_nML2(_ML2),m_nML2_Inh(_INHIB),m_nML2_Inh_Isolated(_INHIB_NOCOUPLING),m_nNeuron(_Connor+_HH+_ML1+_ML2+_INHIB+_INHIB_NOCOUPLING){
//	printf("CoupleVector is %p\n",CoupleVector);
	/*dynamically allocate memory for coupling matrix*/

	m_gProperty=new u8[m_nNeuron];
	memset(m_gProperty,0,sizeof(u8)*m_nNeuron);
	m_sName=new char[40];
	memset(m_sName,0,sizeof(char)*40);
	if(sm_gCouple==NULL){
		sm_gCouple=new double[m_nNeuron];
		memset(sm_gCouple,0,sizeof(double)*m_nNeuron);
	}
	if(sm_gNoise==NULL){
		sm_gNoise=new double[m_nNeuron];
		memset(sm_gNoise,0,sizeof(double)*m_nNeuron);
	}

	m_gNeuron=new neuron*[m_nNeuron];
	int part1=_Connor;
	int part2=_Connor+_HH;
	int part3=_Connor+_HH+_ML1;
	int part4=_Connor+_HH+_ML1+_ML2;
	int part5=_Connor+_HH+_ML1+_ML2+_INHIB;
	int *IsFilled=new int[m_nNeuron];
	memset(IsFilled,0,sizeof(int)*m_nNeuron);
	int index;
	srand((unsigned int)100);
	for(int i=0;i<part1;++i){
		index=rand()%m_nNeuron;
		while(IsFilled[index])
			index=rand()%m_nNeuron;
		IsFilled[index]=1;
		m_gNeuron[index]=new Connor(index,1);
		m_gProperty[index]=NF_CLASS_ONE|NF_RESTING|NF_ISCOUPLED|NF_EXCITATION;
	}
	for(int i=part1;i<part2;++i){
		index=rand()%m_nNeuron;
		while(IsFilled[index])
			index=rand()%m_nNeuron;
		IsFilled[index]=1;
		m_gNeuron[index]=new HH(index,1);
		m_gProperty[i]=NF_CLASS_TWO|NF_RESTING|NF_ISCOUPLED|NF_EXCITATION;
	}
	for(int i=part2;i<part3;++i){
		index=rand()%m_nNeuron;
		while(IsFilled[index])
			index=rand()%m_nNeuron;
		IsFilled[index]=1;
		m_gNeuron[index]=new MorrisLecar(index,1,1);
		m_gProperty[index]=NF_CLASS_ONE|NF_RESTING|NF_ISCOUPLED|NF_EXCITATION;
	}
	for(int i=part3;i<part4;++i){
		index=rand()%m_nNeuron;
		while(IsFilled[index])
			index=rand()%m_nNeuron;
		IsFilled[index]=1;
		m_gNeuron[index]=new MorrisLecar(index,0,1);
		m_gProperty[index]=NF_CLASS_TWO|NF_RESTING|NF_ISCOUPLED|NF_EXCITATION;
	}
	for(int i=part4;i<part5;++i){
		index=rand()%m_nNeuron;
		while(IsFilled[index])
			index=rand()%m_nNeuron;
		IsFilled[index]=1;
		m_gNeuron[index]=new MorrisLecar(index,0,0);
		m_gProperty[index]=(NF_CLASS_TWO|NF_FIRING|NF_INHIBITION)&(~NF_ISCOUPLED);
	}
	for(int i=0;i<m_nNeuron;++i){
		if(!IsFilled[i]){
			m_gNeuron[i]=new MorrisLecar(i,0,0);
			m_gProperty[i]=(NF_CLASS_TWO|NF_FIRING|NF_NOCOUPLING)&(~NF_ISCOUPLED);
		}
	}
	delete[] IsFilled;

}

Network::~Network(){
	/*deallocate memory for coupling matrix*/

	printf("NetworkMatrix is deleted\n");
	delete [] m_gProperty;
	printf("m_gProperty is deleted\n");

	delete []sm_gCouple;
	sm_gCouple=NULL;
	printf("CoupleVector is deleted\n");
	delete []sm_gNoise;
	sm_gNoise=NULL;
	printf("NoiseVector is deleted\n");

	for(int i=0;i<m_nNeuron;++i)
		delete m_gNeuron[i];
	delete []m_gNeuron;
	printf("m_gNeuron is deleted\n");
}


void Network::SetNoiseIntensity(const double _noise_intensity){
	m_noiseintensity=_noise_intensity;
}

//Preprocess 3:Specify neural currents

void Network::ConfigureI(double _ML1_I_Span,double _ML2_I_Span){
	m_ML1_I_Begin=39.7-_ML1_I_Span;
	m_ML1_I_Span=_ML1_I_Span;
	m_ML2_I_Begin=88.1-_ML2_I_Span;
	m_ML2_I_Span=_ML2_I_Span;
}
void Network::RandI(){
	srand((unsigned int)100);
	double ran;
	for(int i=0;i<m_nNeuron;++i){
		ran=Uniform_01();
		if(__IsML1(i)){
			m_gNeuron[i]->SetRandI(ran,m_ML1_I_Begin,m_ML1_I_Span);
		}
		if(__IsML2(i)){
			m_gNeuron[i]->SetRandI(ran,m_ML2_I_Begin,m_ML2_I_Span);
		}				
	}
}
//void Network::RandI(const double _aver_I,const double _deviation){
//	srand((unsigned int)100);
//	double ran1,ran2;
//	for(int i=0;i<m_nNeuron/2;++i){
//		ran1=Uniform_01();
//		ran2=Uniform_01();
//		m_gNeuron[i]->I=_aver_I+sqrt(_deviation)*sqrt(-2*log(ran1))*cos(2*PI*ran2);
//		m_gNeuron[i+m_nNeuron/2]->I=_aver_I+sqrt(_deviation)*sqrt(-2*log(ran1))*sin(2*PI*ran2);
//	}
//	if(m_nNeuron&0x1){
//		ran1=Uniform_01();
//		ran2=Uniform_01();
//		m_gNeuron[m_nNeuron-1]->I=_aver_I+sqrt(_deviation)*sqrt(-2*log(ran1))*cos(2*PI*ran2);
//	}
//}

//Preprocess 4:Set the initial values of variables to sychronization...
void Network::__DoSyncVar(){
	m_t=0.0;

	for(int i=0;i<m_nNeuron;++i)
		m_gNeuron[i]->InitVar();
	srand(100);
}

void Network::__DoAlSyncVar(){
	m_t=0.0;
	srand((unsigned int)100);               //random generator seed is 100
	for(int i=0;i<m_nNeuron;++i)
		m_gNeuron[i]->RandVar();
	srand(100);
}


void Network::SetDt(const double&_dt){
	sm_dt=_dt;
}



void Network::__EulerIterate(){
	for(int i=0;i<m_nNeuron;++i){
		m_gNeuron[i]->euler();
	}
	m_t+=sm_dt;
}

void Network::__RungekuttaIterate(){
	for(int i=0;i<m_nNeuron;++i){
		m_gNeuron[i]->rungekutta();
	}
	m_t+=sm_dt;
}

bool Network::__IsML1(int No){
	if((m_gProperty[No]&NF_CLASS_ONE)&&(m_gProperty[No]&NF_RESTING)){
		if(0==strcmp(m_gNeuron[No]->neuron_model,"MorrisLecar1")){
			return true;
		}
	}
	return false;
}

bool Network::__IsML2(int No){
	if((m_gProperty[No]&NF_CLASS_TWO)&&(m_gProperty[No]&NF_RESTING)){
		if(0==strcmp(m_gNeuron[No]->neuron_model,"MorrisLecar2")){
			return true;
		}
	}
	return false;
}

bool Network::__IsInner(int No){
	int dimension=sqrt(m_nNeuron)+0.5;
	int row=No/dimension;
	int column=No%dimension;
	if(row>=16&&row<=111&&column>=16&&column<=111){
		return true;
	}
	return false;
}

void Network::OutputI(){
	//Create directory
	char pre[50];
	sprintf_s(pre,"F:\\output\\Connor%dHH%dML1_%dML2_%dInh_(%d,%d)",m_nConnor,m_nHH,m_nML1,m_nML2,m_nML2_Inh,m_nML2_Inh_Isolated);
//	_mkdir(pre);

	//Create a file pointer
	FILE *pOutputI;
	char filename[100];
	sprintf_s(filename,"%s_RandI(%.5lf,%.5lf).txt",pre,m_ML1_I_Span,m_ML2_I_Span);
	fopen_s(&pOutputI,filename,"w");
//	fprintf(pOutputI,"NeuronID,ML1_I,ML2_I\n");

	for(int i=0;i<m_nNeuron;++i){
		if(__IsML1(i)){
				fprintf(pOutputI,"%d,%lf,\n",i,m_gNeuron[i]->I);
		}
		if(__IsML2(i)){
				fprintf(pOutputI,"%d,,%lf\n",i,m_gNeuron[i]->I);
		}
	}
}

void Network::OutputNoForOneAndTwo(){
	srand(100);
	FILE *fp;
	char filename[100];
	sprintf_s(filename,"F:\\config\\ML1=%d_ML2=%d_NO.dat",m_nML1,m_nML2);
	fopen_s(&fp,filename,"w");
	int one[5],two[5];
	int indexForOne=0,indexForTwo=0;
	int randNo;
	const int num=5;
	if(m_nML2==0){
		while(indexForOne<num){
			randNo=rand()%m_nNeuron;
			if(__IsInner(randNo)){
				if(__IsML1(randNo)){
					one[indexForOne]=randNo;
					indexForOne++;
				}
			}	
		}
		for(int i=0;i<num;++i){
			fprintf(fp,"%d ",one[i]);
		}
		fprintf(fp,"\n");
	}
	else if(m_nML1==0){
		while(indexForTwo<num){
			randNo=rand()%m_nNeuron;
			if(__IsInner(randNo)){
				if(__IsML2(randNo)){
					two[indexForTwo]=randNo;
					indexForTwo++;
				}
			}		
		}
		for(int i=0;i<num;++i){
			fprintf(fp,"%d ",two[i]);
		}
		fprintf(fp,"\n");
	}
	else{
		while(indexForOne<num&&indexForTwo<num){
			randNo=rand()%m_nNeuron;
			if(__IsInner(randNo)){
				if(__IsML1(randNo)){
					one[indexForOne]=randNo;
					indexForOne++;
				}
				if(__IsML2(randNo)){
					two[indexForTwo]=randNo;
					indexForTwo++;
				}
			}
		}
		while(indexForOne<num){
			randNo=rand()%m_nNeuron;
			if(__IsInner(randNo)){
				if(__IsML1(randNo)){
					one[indexForOne]=randNo;
					indexForOne++;
				}
			}	
		}
		while(indexForTwo<num){
			randNo=rand()%m_nNeuron;
			if(__IsInner(randNo)){
				if(__IsML2(randNo)){
					two[indexForTwo]=randNo;
					indexForTwo++;
				}
			}		
		}
		for(int i=0;i<num;++i){
			fprintf(fp,"%d ",one[i]);
		}
		fprintf(fp,"\n");
		for(int i=0;i<num;++i){
			fprintf(fp,"%d ",two[i]);
		}
		fprintf(fp,"\n");
		
	}
	fclose(fp);
}

void Network::OutputSpikingIndex(){
	//Initialize
	__DoSyncVar();
	//Create directory
	char *direct=new char[100];
	char *specification=new char[40];
	char sSpikingIndex[200];
	__DoMakeDirectory(direct);
	__DoMakeFilename(specification);

	sprintf_s(sSpikingIndex,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_SpikingIndex.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span);	
	/*Create file pointer...*/
	FILE* pOutputSpikingIndex;
	fopen_s(&pOutputSpikingIndex,sSpikingIndex,"w");
	//	fopen_s(&pOutputCouple,coupleFileName,"w");
	int iter_trans=(int)(2000/sm_dt+0.5);
	int iter_end=(int)(5000/sm_dt+0.5);
	//Calculate
	for(int i=1;i<=iter_trans;++i){
		__DoUpdateCouple();
		__UpdateNoise();
		__EulerIterate();
	}
	double *V1=new double[m_nNeuron];
	double *V2=new double[m_nNeuron];
	vector<int> spikingindex;
	
	vector<vector<int> > vecforspikingindex(m_nNeuron);
	for(int i=0;i<m_nNeuron;i++){
		vecforspikingindex[i].push_back(iter_trans+1);
	}
	for(int i=0;i<m_nNeuron;i++){
		V1[i]=m_gNeuron[i]->V;
		V2[i]=m_gNeuron[i]->V;
	}
	for(int i=iter_trans+1;i<=iter_end;++i){
		__DoUpdateCouple();
		__UpdateNoise();
		__EulerIterate();
		
		for(int j=0;j<m_nNeuron;++j){
			if(V2[j]>25&&V2[j]>m_gNeuron[j]->V&&V2[j]>V1[j]){
				vecforspikingindex[j].push_back(i);
			}
			V1[j]=V2[j];
			V2[j]=m_gNeuron[j]->V;
		}
	}
	
	for(int i=0;i<m_nNeuron;i++){
		vecforspikingindex[i].push_back(iter_end);
	}

//	int dimension=sqrt(m_nNeuron)+0.5;
//	int row,column;
	for(int i=0;i<m_nNeuron;++i){
//		row=m_nNeuron/dimension;
//		column=m_nNeuron%dimension;
//		if((row!=0)&&(row!=dimension-1)&&(column!=0)&&(column!=dimension-1)){
			for(auto it=vecforspikingindex[i].begin();it!=vecforspikingindex[i].end();it++){
				fprintf(pOutputSpikingIndex,"%d ",*it);
			}
			fprintf(pOutputSpikingIndex,"\n");
//		}
	}
	fclose(pOutputSpikingIndex);



}

void Network::OutputPhaseAmplitude(){
	//Create directory
	char *direct=new char[100];
	char *specification=new char[40];
	char sSpikingIndex[200],sPhaseAmplitude[200];
	__DoMakeDirectory(direct);
	__DoMakeFilename(specification);

	sprintf_s(sSpikingIndex,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_SpikingIndex.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span);	
	sprintf_s(sPhaseAmplitude,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_PhaseAmplitude.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span);
	ifstream inputSpikingIndex(sSpikingIndex,ios::in);
	string line;
	int iSpikingIndex;
	getline(inputSpikingIndex,line);
	vector<int> maxiList; 
	int maxiListCount;
	istringstream findbeginandend(line);
	while(!findbeginandend.eof()){
		findbeginandend>>iSpikingIndex;
		if(!findbeginandend.fail()){
			maxiList.push_back(iSpikingIndex);
		}
	}
	maxiListCount=maxiList.size();
	int begin=maxiList[0];
	int end=maxiList[maxiListCount-1];
	const int count=end-begin+1;
	int nodecount=0;
	double *cosphase=new double[count];
	double *sinphase=new double[count];
	double *order=new double[count];
	for(int i=0;i<count;++i){
		cosphase[i]=0.0;
		sinphase[i]=0.0;
		order[i]=0.0;
	}

	double deltat,aver_deltat;
	int leftindex,rightindex;
	double phaseangle;
	
	inputSpikingIndex.seekg(0);
	while(!inputSpikingIndex.eof()){
		getline(inputSpikingIndex,line);
		maxiList.clear();
		if(!inputSpikingIndex.fail()){
			nodecount++;
			istringstream instring(line);
			while(!instring.eof()){
				instring>>iSpikingIndex;
				if(!instring.fail()){
					maxiList.push_back(iSpikingIndex);
				}
			}
			maxiListCount=maxiList.size();
			if(maxiListCount<4){
				printf("There are not spikes!\n");
				for(int i=0;i<count;++i){
					cosphase[i]=cosphase[i]+1;
				}
//				return;
				//exit(-1);
			}
			else{
				aver_deltat=double(maxiList[maxiListCount-2]-maxiList[1])/(maxiListCount-3);
				leftindex=maxiList[0];
				rightindex=maxiList[1];
				for(int i=leftindex;i!=rightindex;i++){
					phaseangle=2.0*PI*(aver_deltat+i-rightindex)/aver_deltat;
					cosphase[i-begin]+=cos(phaseangle);
					sinphase[i-begin]+=sin(phaseangle);
				}
 
				for(int m=1;m!=maxiListCount-2;++m){
					leftindex=maxiList[m];
					rightindex=maxiList[m+1];
					deltat=double(rightindex-leftindex);
					for(int i=leftindex;i!=rightindex;++i){
						phaseangle=2.0*PI*(i-leftindex)/deltat;
						cosphase[i-begin]+=cos(phaseangle);
						sinphase[i-begin]+=sin(phaseangle);
					}
				}


                
				leftindex=maxiList[maxiListCount-2];
				rightindex=maxiList[maxiListCount-1];
				for(int i=leftindex;i<=rightindex;++i){
					phaseangle=2.0*PI*(i-leftindex)/aver_deltat;
					cosphase[i-begin]+=cos(phaseangle);
					sinphase[i-begin]+=sin(phaseangle);
				}
			}
            
		}
	}
	for(int i=0;i<count;++i){
		order[i]=sqrt(pow(cosphase[i],2)+pow(sinphase[i],2))/nodecount;
	}

	FILE *pPhaseAmplitude;
	fopen_s(&pPhaseAmplitude,sPhaseAmplitude,"w");
	for(int i=0;i<count;++i){
		fprintf(pPhaseAmplitude,"%lf\n",order[i]);
	}
	fclose(pPhaseAmplitude);
}

void Network::OutputAverISI(){
	//Create directory
	char *direct=new char[100];
	char *specification=new char[40];
	char sSpikingIndex[200],sAverISI[200];
	__DoMakeDirectory(direct);
	__DoMakeFilename(specification);

	FILE *pAverISI;
	sprintf_s(sSpikingIndex,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_SpikingIndex.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span);	
	sprintf_s(sAverISI,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_AverISI.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span);
	fopen_s(&pAverISI,sAverISI,"w");
	ifstream inputSpikingIndex(sSpikingIndex,ios::in);
	string line;
	int iSpikingIndex;
	int iNodeIndex=0;
	double fAverISI;
	vector<int> maxiList; 
	int maxiListCount;
	vector<double> vAverISI;

	
	while(!inputSpikingIndex.eof()){
		getline(inputSpikingIndex,line);
		maxiList.clear();
		if(!inputSpikingIndex.fail()){
			istringstream instring(line);
			while(!instring.eof()){
				instring>>iSpikingIndex;
				if(!instring.fail()){
					maxiList.push_back(iSpikingIndex);
				}
			}
			maxiListCount=maxiList.size();
			if(maxiListCount<4){
				printf("There are not spikes!\n");
				vAverISI.push_back(DBL_MAX);
				continue;
				//exit(-1);
			}
			fAverISI=(double)(maxiList[maxiListCount-2]-maxiList[1])*sm_dt/(maxiListCount-3);
			vAverISI.push_back(fAverISI);      
		}
	}
	for(auto it=vAverISI.begin();it!=vAverISI.end();++it){
		fprintf(pAverISI,"%lf\n",*it);
	}
	fclose(pAverISI);
}

void Network::OutputPopulationFirings(){
	//Create directory
	char *direct=new char[100];
	char *specification=new char[40];
	char sSpikingIndex[200],sPopulationFirings[200];
	__DoMakeDirectory(direct);
	__DoMakeFilename(specification);

	FILE *pPopulationFirings;
	sprintf_s(sSpikingIndex,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_SpikingIndex.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span);	
	sprintf_s(sPopulationFirings,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_PopulationFirings.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span);
	fopen_s(&pPopulationFirings,sPopulationFirings,"w");
	ifstream inputSpikingIndex(sSpikingIndex,ios::in);
	string line;
	int iSpikingIndex;
	int iNodeIndex=0,iNodeRow,iNodeColumn;
	int iDimension=(int)(sqrt(m_nNeuron)+0.5);

	vector<int> maxiList; 
	int maxiListCount;

	while(!inputSpikingIndex.eof()){
		getline(inputSpikingIndex,line);
		maxiList.clear();
		if(!inputSpikingIndex.fail()){
			iNodeRow=iNodeIndex/iDimension;
			iNodeColumn=iNodeIndex%iDimension;
			iNodeIndex++;
			istringstream instring(line);
			while(!instring.eof()){
				instring>>iSpikingIndex;
				if(!instring.fail()){
					maxiList.push_back(iSpikingIndex);
				}
			}
			maxiListCount=maxiList.size();
			if(maxiListCount<4){
				printf("There are not spikes!\n");
				continue;
				//exit(-1);
			}
			double index2t;
			for(int it=1;it!=maxiListCount-1;++it){
				index2t=(double)maxiList[it]*sm_dt;
				fprintf(pPopulationFirings,"%d %d %lf\n",iNodeRow,iNodeColumn,index2t);
			}    
		}
	}
	
	fclose(pPopulationFirings);
}

void Network::OutputPopulationFiringsOnce(){
	//Create directory
	char *direct=new char[100];
	char *specification=new char[40];
	char sSpikingIndex[200],sPopulationFirings[200];
	__DoMakeDirectory(direct);
	__DoMakeFilename(specification);

	FILE *pPopulationFirings;
	sprintf_s(sSpikingIndex,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_SpikingIndex.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span);	
	sprintf_s(sPopulationFirings,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_PopulationFiringsOnce.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span);
	fopen_s(&pPopulationFirings,sPopulationFirings,"w");
	ifstream inputSpikingIndex(sSpikingIndex,ios::in);
	string line;
	int iSpikingIndex;
	int iNodeIndex=0,iNodeRow,iNodeColumn;
	int iDimension=(int)(sqrt(m_nNeuron)+0.5);

	vector<int> maxiList; 
	int maxiListCount;

	while(!inputSpikingIndex.eof()){
		getline(inputSpikingIndex,line);
		maxiList.clear();
		if(!inputSpikingIndex.fail()){
			iNodeRow=iNodeIndex/iDimension;
			iNodeColumn=iNodeIndex%iDimension;
			iNodeIndex++;
			istringstream instring(line);
			while(!instring.eof()){
				instring>>iSpikingIndex;
				if(!instring.fail()){
					maxiList.push_back(iSpikingIndex);
				}
			}
			maxiListCount=maxiList.size();
			if(maxiListCount<4){
				printf("There are not spikes!\n");
				continue;
				//exit(-1);
			}
			double index2t=(double)maxiList[1]*sm_dt;
			fprintf(pPopulationFirings,"%d %d %lf\n",iNodeRow,iNodeColumn,index2t);
	 
		}
	}	
	fclose(pPopulationFirings);
}

void Network::OutputNoise(){
	FILE *pNoise;
	fopen_s(&pNoise,"Noise.dat","w");
	for(int i=0;i<100000;++i){
		__UpdateNoise();
		fprintf(pNoise,"%lf\n",sm_gNoise[0]);
	}
	fclose(pNoise);
}











/*Update noise*/
void Network::__UpdateNoise(){
	if(abs(m_noiseintensity)<0.0000001)
		return;
	double ran1,ran2;
	double log_temp,angle_temp;
	for(int i=0;i<m_nNeuron;i+=2){  //m_nNeuron must be even.
		ran1=Uniform_01();
		ran2=Uniform_01();
		log_temp=log(ran1);
		angle_temp=2.0*PI*ran2;
		//(m_noiseintensity/sqrt(sm_dt))*//*sqrt(2.0*m_noiseintensity)
		sm_gNoise[i]=sqrt(-4.0*m_noiseintensity*log_temp)*cos(angle_temp);
		//(m_noiseintensity/sqrt(sm_dt))*//*sqrt(2.0*m_noiseintensity)
		sm_gNoise[i+1]=sqrt(-4.0*m_noiseintensity*log_temp)*sin(angle_temp);
	}
	//if(m_nNeuron&0x1){
	//	ran1=Uniform_01();
	//	ran2=Uniform_01();
	//	NoiseVector[m_nNeuron-1]=sqrt(2.0*m_noiseintensity)*sqrt(-2.0*log(ran1))*cos(2.0*PI*ran2);
	//}
}

double Network::SyncFactor(){
	double *sum_VikSquare=new double[m_nNeuron];
	double *sum_Vik=new double[m_nNeuron];
	double sum_averVkSquare=0.0,sum_averVk=0.0;

	double sumVecField=0.0,averVecField=0.0;
	double numerator=0.0,denominator=0.0,R=0.0;
	const int iter_end=40001;
	const int eval_begin=20000;
	const int nTime=iter_end-eval_begin-1;
	memset(sum_VikSquare,0,sizeof(double)*m_nNeuron);
	memset(sum_Vik,0,sizeof(double)*m_nNeuron);
	__DoAlSyncVar();
	//Calculate
	for(int iter=1;iter<iter_end;++iter){
	 	__DoUpdateCouple();
	//	__UpdateNoise();

		__EulerIterate();
		if(iter>eval_begin){
			sumVecField=0.0;
			for(int index=0;index<m_nNeuron;++index){
				sumVecField+=m_gNeuron[index]->V;
				sum_Vik[index]+=m_gNeuron[index]->V;
				sum_VikSquare[index]+=(m_gNeuron[index]->V)*(m_gNeuron[index]->V);
			}
			averVecField=sumVecField/m_nNeuron;
			sum_averVk+=(averVecField);
			sum_averVkSquare+=(averVecField*averVecField);
		}
	}

	numerator=((sum_averVkSquare/nTime)-(sum_averVk/nTime)*(sum_averVk/nTime));
	for(int index=0;index<m_nNeuron;++index){
			denominator+=((sum_VikSquare[index]/nTime)-(sum_Vik[index]/nTime)*(sum_Vik[index]/nTime));
	}
	denominator/=m_nNeuron;
	R=numerator/denominator;

	delete []sum_Vik;
	delete []sum_VikSquare;

	return R;

}

//double Network::OrderParameter(const double _eval_begin,const double _iter_end){
//	//Initialize
//	__DoAlSyncVar();
//	SetDt(0.01);
//	const int eval_begin=(int)(_eval_begin/sm_dt);
//	const int iter_end=(int)(_iter_end/sm_dt);
//	int nTime=0;
//
//	//calculate
//	for(int iter=1;iter!=eval_begin;++iter){
//		UpdateChemical();
//		__RungekuttaIterate();
//		UpdateSpike();
//	}
//	int k;
//	//memset(spike_num,1,sizeof(int)*m_nNeuron);
//	for(int index=0;index!=m_nNeuron;++index){
//		if(spike_time[index][0]>0) spike_num[index]=1;
//		else spike_num[index]=0;
//	}
//	for(int iter=eval_begin;iter!=iter_end;++iter){
//		memset(sm_gCouple,0,sizeof(double)*m_nNeuron);
//		for(int i=0;i<m_nNeuron;++i){
//			for(int j=0;j<m_nNeuron;++j){
//				if(NetworkMatrix[i][j]){
//					k=spike_num[j]-1;
//					CoupleVector[i]+=NormalA*(exp(-(m_t-spike_time[j][k])/tau1)-exp(-(m_t-spike_time[j][k])/tau2));
//				}
//			}
//			CoupleVector[i]=CoupleVector[i]*gc_c*(Vsyn_exc-m_gNeuron[i]->V);
//	//		printf("CoupleVector[%d] is %lf\n",i,CoupleVector[i]);
//		}
//		__RungekuttaIterate();
//		for(int index=0;index<m_nNeuron;++index){
//			if(V_pre[index]>0.5&&V_pre[index]>m_gNeuron[index]->V&&V_pre[index]>V_pre1[index]){
//				if(spike_num[index]!=m_nNeuron){
//					spike_time[index][spike_num[index]]=m_t;
//					spike_num[index]=spike_num[index]+1;
//				}
//				else{
//					memmove(spike_time[index],spike_time[index]+1,sizeof(double)*(m_nNeuron-1));
//					spike_time[index][m_nNeuron-1]=m_t;
//				}
//			}
//			V_pre1[index]=V_pre[index];
//			V_pre[index]=m_gNeuron[index]->V;
//		}
//	}
///*	int firing_end;
//	for(int index=0;index!=m_nNeuron;++index){
//		firing_end=spike_num[index];
//		printf("neuron %d %s spike time: ",index,m_gNeuron[index]->neuron_model);
//		for(int firing_num=0;firing_num!=spike_num[index];++firing_num){
//			printf(" %lf",spike_time[index][firing_num]);
//		}
//		printf("\n");
//	}*/
//	double feval_end=spike_time[0][spike_num[0]-1];
//	for(int index=1;index!=m_nNeuron;++index){
//		if(feval_end>spike_time[index][spike_num[index]-1]) feval_end=spike_time[index][spike_num[index]-1];
//	}
////	printf("feval_end is %lf\n",feval_end);
//	int ieval_end=(int)(feval_end/sm_dt);
//	double R=0.0;
//	double re=0.0,im=0.0;
//	double phi;
//	double *pPeriod=new double[m_nNeuron];
//	double *T_1k=new double[m_nNeuron];
//	double *T_2k=new double[m_nNeuron];
//	double eval_t;
//	for(int index=0;index!=m_nNeuron;++index){
//		T_1k[index]=spike_time[index][0];
//		T_2k[index]=spike_time[index][1];
//		pPeriod[index]=T_2k[index]-T_1k[index];
//	}
//		
//	memset(spike_num,0,sizeof(int)*m_nNeuron);
//	for(int iter=eval_begin;iter<ieval_end;iter++){
//		eval_t=sm_dt*iter;
//		for(int index=0;index!=m_nNeuron;++index){
//			if(T_2k[index]<=eval_t){
//				k=spike_num[index]+1;
//				spike_num[index]=k;
//				T_1k[index]=spike_time[index][k];
//				T_2k[index]=spike_time[index][k+1];
//				pPeriod[index]=T_2k[index]-T_1k[index];
//			}
//		}
//		re=0.0,im=0.0;
//		for(int index=0;index<m_nNeuron;++index){
//			phi=2.0*PI*(eval_t-T_1k[index])/pPeriod[index];
//			re+=cos(phi);
//			im+=sin(phi);
//		}
//	    R+=sqrt(pow(re,2.0)+pow(im,2.0));
//		++nTime;
//	//	printf("%d\n",iter);
//	}
////	printf("nTime is %d\n",nTime);
//	R/=m_nNeuron*nTime;
//	
//	delete []pPeriod;
//	delete []T_1k;
//	delete []T_2k;
//
//	return R;
//}
//process
void Network::OutputTimeSeries(){
	//Initialize
	__DoSyncVar();
	//Create directory
	char *direct=new char[100];
	char *specification=new char[40];
	char timeSeriesFileName[200];

	__DoMakeDirectory(direct);
	__DoMakeFilename(specification);

	sprintf_s(timeSeriesFileName,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_TimeSeries.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span);
//	sprintf_s(coupleFileName,"%s\\%s_noise=%.5lf_Couple.dat",direct,specification,m_noiseintensity);
	/*Create file pointer...*/
	FILE* pOutputTimeSeries;
	fopen_s(&pOutputTimeSeries,timeSeriesFileName,"w");
//	fopen_s(&pOutputCouple,coupleFileName,"w");
	int iter_trans=(int)(1000/sm_dt+0.5);
	int iter_end=(int)(3000/sm_dt+0.5);
	//Calculate

	char inFilename[100];
	sprintf_s(inFilename,"F:\\config\\ML1=%d_ML2=%d_NO.dat",m_nML1,m_nML2);
	std::ifstream No(inFilename,std::ios::in);
	std::vector<int> vec;
	int inputNo;
	while(!No.eof()){
		No>>inputNo;
		if(No.fail())
			break;
		vec.push_back(inputNo);
	}
	No.close();
	int n=vec.size();
	for(int i=1;i<=iter_trans;++i){
		__DoUpdateCouple();
		__UpdateNoise();
		__EulerIterate();
	}
	for(int i=iter_trans+1;i<=iter_end;++i){
		__DoUpdateCouple();
		__UpdateNoise();
		__EulerIterate();
		fprintf(pOutputTimeSeries,"%lf ",m_t);
		for(int j=0;j<n;++j){
			fprintf(pOutputTimeSeries,"%lf ",m_gNeuron[vec[j]]->V);
		}
		fprintf(pOutputTimeSeries,"\n");
		//printf("%d\n",i);
	}
	fclose(pOutputTimeSeries);
}

void Network::OutputCouple(){
	//Initialize
	__DoAlSyncVar();
	//Create directory
	char *direct=new char[100];
	char *specification=new char[40];
	char coupleCurrentFileName[200];

	__DoMakeDirectory(direct);
	__DoMakeFilename(specification);

	sprintf_s(coupleCurrentFileName,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_CoupleCurrent.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span);
//	sprintf_s(coupleFileName,"%s\\%s_noise=%.5lf_Couple.dat",direct,specification,m_noiseintensity);
	/*Create file pointer...*/
	FILE* pOutputCoupleCurrent;
	fopen_s(&pOutputCoupleCurrent,coupleCurrentFileName,"w");
//	fopen_s(&pOutputCouple,coupleFileName,"w");
	int iter_trans=(int)(2000/sm_dt+0.5);
	int iter_end=(int)(3000/sm_dt+0.5);
	//Calculate

	char inFilename[100];
	sprintf_s(inFilename,"F:\\config\\ML1=%d_ML2=%d_NO.dat",m_nML1,m_nML2);
	std::ifstream No(inFilename,std::ios::in);
	std::vector<int> vec;
	int inputNo;
	while(!No.eof()){
		No>>inputNo;
		if(No.fail())
			break;
		vec.push_back(inputNo);
	}
	No.close();
	int n=vec.size();
	for(int i=1;i<=iter_trans;++i){
		__DoUpdateCouple();
		__UpdateNoise();
		__EulerIterate();
	}
	for(int i=iter_trans+1;i<=iter_end;++i){
		__DoUpdateCouple();
		__UpdateNoise();
		__EulerIterate();
		fprintf(pOutputCoupleCurrent,"%lf ",m_t);
		for(int j=0;j<n;++j){
			fprintf(pOutputCoupleCurrent,"%lf %lf ",m_gNeuron[vec[j]]->V,sm_gCouple[vec[j]]);
		}
		fprintf(pOutputCoupleCurrent,"\n");
		//printf("%d\n",i);
	}
	fclose(pOutputCoupleCurrent);
}
void Network::SpiralWave(){
	/*Initialize...*/
	__DoSyncVar();
	/*Create directory...*/

	char *direct=new char[100];
	char *specification=new char[40];
	char filename[200];
	__DoMakeDirectory(direct);
	__DoMakeFilename(specification);

//	int length=strlen(specification);
	int iColumn=int(sqrt(m_nNeuron)+0.5);
	/*Create file pointer...*/
	FILE* pOutputSpiralWave;

	int iDimension=int(sqrt((double)m_nNeuron)+0.5);
	for(int i=0;i<30000;++i){
		if(i>10000&&i%200==0){
			sprintf_s(filename,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_t=%.5lf.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span,m_t);
			//sprintf_s(filename,"%s\\%s_noise=%.5lf_t=%.5lf.dat",direct,specification,m_noiseintensity,m_t);
			fopen_s(&pOutputSpiralWave,filename,"w");
	
			for(int j=0;j<m_nNeuron;++j){
				fprintf(pOutputSpiralWave,"%lf ",m_gNeuron[j]->V);
				if(j%iDimension==iColumn-1) fprintf(pOutputSpiralWave,"\n");
			}
			fprintf(pOutputSpiralWave,"\n");
			fclose(pOutputSpiralWave);
		}
	
		//Calculate...
		__DoUpdateCouple();
     	__UpdateNoise();
		__EulerIterate();
	}
	delete []direct;
	delete []specification;
	
}

void Network::SpiralWave(int _begin,int _end,int _dt){
	/*Initialize...*/
	__DoSyncVar();
	/*Create directory...*/

	char *direct=new char[100];
	char *specification=new char[40];
	char filename[200];
	__DoMakeDirectory(direct);
	__DoMakeFilename(specification);

//	int length=strlen(specification);
	int iColumn=int(sqrt(m_nNeuron)+0.5);
	/*Create file pointer...*/
	FILE* pOutputSpiralWave;

	int iDimension=int(sqrt((double)m_nNeuron)+0.5);
	int iBegin=int(_begin/sm_dt+0.5);
	int iEnd=int(_end/sm_dt+0.5);
	int iDt=int(_dt/sm_dt+0.5);
	for(int i=0;i<=iEnd;++i){
		if(i>=iBegin&&i%iDt==0){
			sprintf_s(filename,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_t=%.5lf.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span,m_t);
			//sprintf_s(filename,"%s\\%s_noise=%.5lf_t=%.5lf.dat",direct,specification,m_noiseintensity,m_t);
			fopen_s(&pOutputSpiralWave,filename,"w");
	
			for(int j=0;j<m_nNeuron;++j){
				fprintf(pOutputSpiralWave,"%lf ",m_gNeuron[j]->V);
				if(j%iDimension==iColumn-1) fprintf(pOutputSpiralWave,"\n");
			}
			fprintf(pOutputSpiralWave,"\n");
			fclose(pOutputSpiralWave);
		}
	
		//Calculate...
		__DoUpdateCouple();
     	__UpdateNoise();
		__EulerIterate();
	}
	delete []direct;
	delete []specification;
	
}

void Network::SpiralWaveAndTimeSeries(){
	/*Initialize...*/
	__DoSyncVar();
	/*Create directory...*/

	char *direct=new char[100];
	char *specification=new char[40];
	char filename[200];
	__DoMakeDirectory(direct);
	__DoMakeFilename(specification);
	/*Create file pointer...*/
	FILE* pOutputSpiralWave,*pOutputTimeSeries;
	sprintf_s(filename,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_TimeSeries.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span);

	fopen_s(&pOutputTimeSeries,filename,"w");
	int iDimension=int(sqrt((double)m_nNeuron)+0.5);
	int iter_trans=(int)(1000/sm_dt+0.5);
	int iter_end=(int)(3000/sm_dt+0.5);
	//Calculate

	char inFilename[100];
	sprintf_s(inFilename,"F:\\config\\ML1=%d_ML2=%d_NO.dat",m_nML1,m_nML2);
	std::ifstream No(inFilename,std::ios::in);
	std::vector<int> vec;
	int inputNo;
	while(!No.eof()){
		No>>inputNo;
		if(No.fail())
			break;
		vec.push_back(inputNo);
	}
	No.close();
	int n=vec.size();
	for(int i=1;i<=iter_trans;++i){
		__DoUpdateCouple();
		__UpdateNoise();
		__EulerIterate();
	}

	for(int i=iter_trans+1;i<iter_end;++i){
		__DoUpdateCouple();
		__UpdateNoise();
		__EulerIterate();
		fprintf(pOutputTimeSeries,"%lf ",m_t);
		for(int j=0;j<n;++j){
			fprintf(pOutputTimeSeries,"%lf ",m_gNeuron[vec[j]]->V);
		}
		fprintf(pOutputTimeSeries,"\n");
		if(i%200==0){
			sprintf_s(filename,"%s\\%s_noise=%.5lf_RandI(%.5lf,%.5lf)_t=%.5lf.dat",direct,specification,m_noiseintensity,m_ML1_I_Span,m_ML2_I_Span,m_t);
		    fopen_s(&pOutputSpiralWave,filename,"w");	

			for(int j=0;j<m_nNeuron;++j){
				fprintf(pOutputSpiralWave,"%lf ",m_gNeuron[j]->V);
				if(j%iDimension==iDimension-1) fprintf(pOutputSpiralWave,"\n");
			}
			fprintf(pOutputSpiralWave,"\n");
			fclose(pOutputSpiralWave);
		}
	}
	
	fclose(pOutputTimeSeries);
	delete []direct;
	delete []specification;
}

void Network::MaximumPotential(){
	/*Initialize...*/
	__DoSyncVar();
	/*Create directory...*/

	char *direct=new char[100];
	char *specification=new char[40];
	char filename[140];
	__DoMakeDirectory(direct);
	__DoMakeFilename(specification);
	int iColumn=int(sqrt(m_nNeuron)+0.5);
	/*Create file pointer...*/
	FILE* pOutputSpiralWave;
	double *aMaximumV=new double[m_nNeuron];
	double *aMinimumV=new double[m_nNeuron];
	for(int i=0;i<m_nNeuron;++i){
		aMaximumV[i]=-DBL_MAX;
	}
	for(int i=0;i<m_nNeuron;++i){
		aMinimumV[i]=DBL_MAX;
	}
	int iDimension=int(sqrt((double)m_nNeuron)+0.5);
	for(int i=0;i<30000;++i){
		if(i>20000&&i%20==0){
			for(int j=0;j<m_nNeuron;++j){
				if(aMaximumV[j]<m_gNeuron[j]->V)
					aMaximumV[j]=m_gNeuron[j]->V;
				if(aMinimumV[j]>m_gNeuron[j]->V)
					aMinimumV[j]=m_gNeuron[j]->V;
			}
		}
	
		//Calculate...
		__DoUpdateCouple();
     	__UpdateNoise();
		__EulerIterate();
	}
	sprintf_s(filename,"%s\\%s_noise=%.5lf_MaximumV.dat",direct,specification,m_noiseintensity);
	fopen_s(&pOutputSpiralWave,filename,"w");
	for(int j=0;j<m_nNeuron;++j){
		fprintf(pOutputSpiralWave,"%lf ",aMaximumV[j]);
		if(j%iDimension==iColumn-1) fprintf(pOutputSpiralWave,"\n");
	}
	fprintf(pOutputSpiralWave,"\n");
	fclose(pOutputSpiralWave);

	sprintf_s(filename,"%s\\%s_noise=%.5lf_MinimumV.dat",direct,specification,m_noiseintensity);
	fopen_s(&pOutputSpiralWave,filename,"w");
	for(int j=0;j<m_nNeuron;++j){
		fprintf(pOutputSpiralWave,"%lf ",aMinimumV[j]);
		if(j%iDimension==iColumn-1) fprintf(pOutputSpiralWave,"\n");
	}
	fprintf(pOutputSpiralWave,"\n");
	fclose(pOutputSpiralWave);

	delete []direct;
	delete []specification;
	delete []aMaximumV;
}
//void Network::Tau1Tau2SyncFactor(const double _tau_max,const double _interval,const double _g){
//	//Create directory...
//	char dir[40];
//	sprintf_s(dir,"Connor%dHH%dML1_%dML2_%d",Connor_num,HH_num,ML1_num,ML2_num);
//	_mkdir(dir);
//	//Create file pointer...
//	FILE *pOutputSyncFactor;
//	char Outfilename[60];
//	sprintf_s(Outfilename,".\\%s\\gc_c=%.4lf_tau1_tau2_sync_factor.dat",dir,_g);
//	fopen_s(&pOutputSyncFactor,Outfilename,"w");
//
//	double *sum_VikSquare=new double[m_nNeuron];
//	double *sum_Vik=new double[m_nNeuron];
//	double sum_averVkSquare=0.0,sum_averVk=0.0;
//
//	double sumVecField,averVecField;
//	double numerator,denominator,R;
//	double averVecFieldMax,averVecFieldMin,amplitude;
//	const int iter_end=40001;
//	const int eval_begin=20000;
//	const int nTime=iter_end-eval_begin;
//	gc_c=_g;
//	int iternum=(int)(_tau_max/_interval);
//	int tau1_beg=(int)(2/_interval);
//	int tau2_beg=(int)(1/_interval);
//	for(int taui=tau1_beg;taui<iternum;++taui){
//		for(int tauj=tau2_beg;tauj<taui;++tauj){
//			//assign parameter
//			tau1=_interval*taui;
//			tau2=_interval*tauj;
//
//			numerator=0.0;
//			denominator=0.0;
//			sumVecField=0.0;
//			sum_averVkSquare=0.0;
//			sum_averVk=0.0;
//			memset(sum_VikSquare,0,sizeof(double)*m_nNeuron);
//			memset(sum_Vik,0,sizeof(double)*m_nNeuron);
//			__DoAlSyncVar();
//
//			for(int iter=1;iter<eval_begin;++iter){
//	 			UpdateChemical();
//				
//				__RungekuttaIterate();
//				UpdateSpike();
//
//			}
//			
//			sumVecField=0.0;
//			for(int index=0;index!=m_nNeuron;++index)
//				sumVecField+=m_gNeuron[index]->V;
//			averVecField=sumVecField/m_nNeuron;
//			averVecFieldMax=averVecFieldMin=averVecField;
//			for(int iter=eval_begin;iter!=iter_end;++iter){
//				UpdateChemical();
//				//UpdateElectrical();
//				__RungekuttaIterate();
//	    		UpdateSpike();
//				sumVecField=0.0;
//				for(int index=0;index<m_nNeuron;++index){
//					sumVecField+=m_gNeuron[index]->V;
//					sum_Vik[index]+=m_gNeuron[index]->V;
//					sum_VikSquare[index]+=(m_gNeuron[index]->V)*(m_gNeuron[index]->V);
//				}
//				averVecField=sumVecField/m_nNeuron;
//				sum_averVk+=(averVecField);
//				sum_averVkSquare+=(averVecField*averVecField);
//				if(averVecFieldMax<averVecField) averVecFieldMax=averVecField;
//				if(averVecFieldMin>averVecField) averVecFieldMin=averVecField;
//			}
//	
//			numerator=((sum_averVkSquare/nTime)-(sum_averVk/nTime)*(sum_averVk/nTime));
//			for(int index=0;index<m_nNeuron;++index){
//					denominator+=((sum_VikSquare[index]/nTime)-(sum_Vik[index]/nTime)*(sum_Vik[index]/nTime));
//			}
//			denominator/=m_nNeuron;
//			R=numerator/denominator;
//			amplitude=averVecFieldMax-averVecFieldMin;
//			//Output
//			fprintf(pOutputSyncFactor,"%lf %lf %lf %lf\n",tau1,tau2,R,amplitude);
//	
//			printf("tau1=%lf tau2=%lf gc_c=%lf R=%lf amplitude=%lf is completed!\n",taui,tauj,gc_c,R,amplitude);
//
//			}
//	}
//	fclose(pOutputSyncFactor);
//
//
//	delete []sum_Vik;
//	delete []sum_VikSquare;
//
//}

//void Network::Tau1Tau2OrderParameter(const double _tau_max,const double _interval,const double _g){
//	//Create directory...
//	char dir[40];
//	sprintf_s(dir,"Connor%dHH%dML1_%dML2_%d",Connor_num,HH_num,ML1_num,ML2_num);
//	_mkdir(dir);
//	//Create file pointer...
//	FILE *pOutputOrderParameter;
//	char Outfilename[80];
//	gc_c=_g;
//	int iternum=(int)(_tau_max/_interval);
//	int tau1_beg=(int)(2/_interval);
//	int tau2_beg=(int)(1/_interval);
//	for(int taui=tau1_beg;taui<iternum;++taui){
//		for(int tauj=tau2_beg;tauj<taui;++tauj){
//			//assign parameter
//			tau1=_interval*taui;
//			tau2=_interval*tauj;
//	
//			sprintf_s(Outfilename,".\\%s\\tau1=%.4lf_tau2=%.4lf_g=%.4lf_order_parameter.dat",dir,tau1,tau2,gc_c);
//			fopen_s(&pOutputOrderParameter,Outfilename,"w");
//			fprintf(pOutputOrderParameter,"tau1=%.4lf_tau2=%.4lf_g=%.4lf\n",tau1,tau2,gc_c);
//			double R;
//			int eval_begin;
//			int eval_end;
//			for(int eval=1;eval<40;++eval){
//				eval_begin=200*eval;
//				eval_end=eval_begin+200;
//				R=OrderParameter(eval_begin,eval_end);
//				fprintf(pOutputOrderParameter,"%lf\n",R);
//			}
//			fclose(pOutputOrderParameter);
//		}
//	}
//}


SparseNetwork::SparseNetwork():Network(){
	m_gCoupleIdx=new int*[m_nNeuron];
	m_gCoupleNum=new int[m_nNeuron];
}

SparseNetwork::SparseNetwork(int _Connor,int _HH,int _ML1,int _ML2,int _Inh,int _InhNoCoupling):Network(_Connor,_HH,_ML1,_ML2,_Inh,_InhNoCoupling){
	m_gCoupleIdx=new int*[m_nNeuron];
	m_gCoupleNum=new int[m_nNeuron];
}

SparseNetwork::~SparseNetwork(){
	delete []m_gCoupleIdx;
	delete []m_gCoupleNum;
}

void SparseNetwork::BuildNearestNeighbor(int _adjacent){
	//Generate the nearest neighbor coupling network
	if(_adjacent<1) fprintf(stderr,"Invalid Input!\n");
	m_gCoupleIdx[0]=new int[m_nNeuron*2*_adjacent];
	m_gCoupleNum=new int[m_nNeuron];
	for(int i=0;i<m_nNeuron;++i){
		m_gCoupleNum[i]=2*_adjacent;
	}
	for(int i=1;i<m_nNeuron;++i){
		m_gCoupleIdx[i]=m_gCoupleIdx[i-1]+2*_adjacent;
	}
	int flag;
	int icr;
	for(int i=0;i<m_nNeuron;++i){
		icr=0;
		for(int j=1;j<_adjacent+1;++j){
			flag=i-j;
			if(flag<0)
				flag=flag+m_nNeuron;
			m_gCoupleIdx[i][icr]=flag;
			icr++;

			flag=i+j;
			if(flag>m_nNeuron)
				flag=flag-m_nNeuron;
			m_gCoupleIdx[i][icr]=flag;
			icr++;
		}
	}
}


void SparseNetwork::BuildSquareSmallWorld(const double _rewiring){
	sprintf_s(m_sName,40,"SmallWorld_%.2f",_rewiring);
	int nDimension=sqrt(m_nNeuron)+0.5;
	bool **RelationMatrix=new bool*[m_nNeuron];
	RelationMatrix[0]=new bool[m_nNeuron*m_nNeuron];
	for(int i=1;i<m_nNeuron;++i){
		RelationMatrix[i]=RelationMatrix[i-1]+m_nNeuron;
	}
	memset(RelationMatrix[0],0,sizeof(bool)*m_nNeuron*m_nNeuron);
	int iRow,iColumn;
	for(int i=0;i<m_nNeuron;++i){
		iRow=i/nDimension;
		iColumn=i%nDimension;
		//up
		if(iRow>0){
			RelationMatrix[i][(iRow-1)*nDimension+iColumn]=true;
			RelationMatrix[(iRow-1)*nDimension+iColumn][i]=true;
		}

		//left
		if(iColumn>0){
			RelationMatrix[i][iRow*nDimension+iColumn-1]=true;
			RelationMatrix[iRow*nDimension+iColumn-1][i]=true;
		}

		//right
		if(iColumn<nDimension-1){
			RelationMatrix[i][iRow*nDimension+iColumn+1]=true;
			RelationMatrix[iRow*nDimension+iColumn+1][i]=true;
		}

		//down
		if(iRow<nDimension-1){
			RelationMatrix[i][(iRow+1)*nDimension+iColumn]=true;
			RelationMatrix[(iRow+1)*nDimension+iColumn][i]=true;
		}
	}

	srand((unsigned int)100);
	int iElement;
	for(int i=0;i<m_nNeuron;++i){
		for(int j=i+1;j<m_nNeuron;++j){
			if(RelationMatrix[i][j]&&(rand()/(double)RAND_MAX)<_rewiring){
				RelationMatrix[i][j]=false;
				RelationMatrix[j][i]=false;
				while((iElement=rand()%m_nNeuron)==i);
				RelationMatrix[i][iElement]=true;
				RelationMatrix[iElement][i]=true;
			}
		}
	}
	
	vector<int> vec_int;
	for(int i=0;i<m_nNeuron;++i){
		for(int j=0;j<m_nNeuron;++j){
			if(RelationMatrix[i][j]){
				vec_int.push_back(j);
			}
		}
		size_t nCount;
		nCount=m_gCoupleNum[i]=vec_int.size();
		
		m_gCoupleIdx[i]=new int[nCount];

		for(int idx=0;idx<nCount;++idx){
			m_gCoupleIdx[i][idx]=vec_int[idx];
		}
		vec_int.clear();
	}
	delete []RelationMatrix[0];
	delete []RelationMatrix;

}

void SparseNetwork::BuildSquareSparser(const double _sparse){
	sprintf_s(m_sName,40,"Sparse_%.2f",_sparse);
	int nDimension=sqrt(m_nNeuron)+0.5;
	bool **RelationMatrix=new bool*[m_nNeuron];
	RelationMatrix[0]=new bool[m_nNeuron*m_nNeuron];
	for(int i=1;i<m_nNeuron;++i){
		RelationMatrix[i]=RelationMatrix[i-1]+m_nNeuron;
	}
	memset(RelationMatrix[0],0,sizeof(bool)*m_nNeuron*m_nNeuron);
	int iRow,iColumn;
	for(int i=0;i<m_nNeuron;++i){
		iRow=i/nDimension;
		iColumn=i%nDimension;
		//up
		if(iRow>0){
			RelationMatrix[i][(iRow-1)*nDimension+iColumn]=true;
			RelationMatrix[(iRow-1)*nDimension+iColumn][i]=true;
		}

		//left
		if(iColumn>0){
			RelationMatrix[i][iRow*nDimension+iColumn-1]=true;
			RelationMatrix[iRow*nDimension+iColumn-1][i]=true;
		}

		//right
		if(iColumn<nDimension-1){
			RelationMatrix[i][iRow*nDimension+iColumn+1]=true;
			RelationMatrix[iRow*nDimension+iColumn+1][i]=true;
		}

		//down
		if(iRow<nDimension-1){
			RelationMatrix[i][(iRow+1)*nDimension+iColumn]=true;
			RelationMatrix[(iRow+1)*nDimension+iColumn][i]=true;
		}
	}

	srand((unsigned int)100);

	for(int i=0;i<m_nNeuron;++i){
		for(int j=i+1;j<m_nNeuron;++j){
			if(RelationMatrix[i][j]&&(rand()/(double)RAND_MAX)<_sparse){
				RelationMatrix[i][j]=false;
				RelationMatrix[j][i]=false;
			}
		}
	}
	
	vector<int> vec_int;
	for(int i=0;i<m_nNeuron;++i){
		for(int j=0;j<m_nNeuron;++j){
			if(RelationMatrix[i][j]){
				vec_int.push_back(j);
			}
		}
		int nCount;
		nCount=m_gCoupleNum[i]=vec_int.size();
		
		m_gCoupleIdx[i]=new int[nCount];

		for(int idx=0;idx<nCount;++idx){
			m_gCoupleIdx[i][idx]=vec_int[idx];
		}
		vec_int.clear();
	}
	delete []RelationMatrix[0];
	delete []RelationMatrix;
}

void SparseNetwork::BuildSquare(){
	sprintf_s(m_sName,40,"Square");
	int iRow,iColumn;
	int iDimension=sqrt(m_nNeuron)+0.5;
	int aTempIdx[4];
	int nTotal;

	for(int i=0;i<m_nNeuron;++i){
		iRow=i/iDimension;
		iColumn=i%iDimension;
        //up
		nTotal=0;
		if((iRow-1)>=0){
			aTempIdx[nTotal]=(iRow-1)*iDimension+iColumn;
			nTotal++;
		}

		//left
		if((iColumn-1)>=0){
			aTempIdx[nTotal]=iRow*iDimension+iColumn-1;
			nTotal++;
		}

		//down
		if((iRow+1)<iDimension){
			aTempIdx[nTotal]=(iRow+1)*iDimension+iColumn;
			nTotal++;
		}

		//right
		if((iColumn+1)<iDimension){
			aTempIdx[nTotal]=iRow*iDimension+iColumn+1;
			nTotal++;
		}
		
		m_gCoupleIdx[i]=new int[nTotal];

		for(int j=0;j<nTotal;++j){
			m_gCoupleIdx[i][j]=aTempIdx[j];
		}

		m_gCoupleNum[i]=nTotal;
	}
}

void SparseNetwork::OutputCoupleIdx(){
	FILE *pSquare;
	fopen_s(&pSquare,"CoupleIdx.dat","w");
	int iDimension=(int)(sqrt(m_nNeuron)+0.5);
	int nCount;
	for(int i=0;i<m_nNeuron;++i){
//		fprintf(pSquare,"neuron(%d)(Exc %d,Inh %d) ",i,aExcitation[i],aInhibition[i]);
		fprintf(pSquare,"neuron(%d) ",i);
		if(m_gProperty[i]&NF_ISCOUPLED){
			nCount=m_gCoupleNum[i];//aExcitation[i]+aInhibition[i];
			for(int j=0;j<nCount;++j){
				fprintf(pSquare,"%d ",m_gCoupleIdx[i][j]);
			}
		}
		fprintf(pSquare,"\n");
	}
	fclose(pSquare);
}