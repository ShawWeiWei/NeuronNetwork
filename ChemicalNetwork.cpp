#include "stdafx.h"

#include "ChemicalNetwork.h"


ChemicalSNetwork::ChemicalSNetwork():SparseNetwork(){
	aExcitation=new int[m_nNeuron];
	memset(aExcitation,0,sizeof(int)*m_nNeuron);
	aInhibition=new int[m_nNeuron];
	memset(aInhibition,0,sizeof(int)*m_nNeuron);
}

ChemicalSNetwork::ChemicalSNetwork(int _Connor,int _HH,int _ML1,int _ML2,int _Inh,int _InhNoCoupling):SparseNetwork(_Connor,_HH,_ML1,_ML2,_Inh,_InhNoCoupling){
	aExcitation=new int[m_nNeuron];
	memset(aExcitation,0,sizeof(int)*m_nNeuron);
	aInhibition=new int[m_nNeuron];
	memset(aInhibition,0,sizeof(int)*m_nNeuron);
}

ChemicalSNetwork::~ChemicalSNetwork(){
	delete []aExcitation;
	delete []aInhibition;
}

/*up*/
/*if((iRow-1)>=0){
	iTempIndex=(iRow-1)*iDimension+iColumn;
	if(m_gProperty[iTempIndex]&NF_EXCITATION){
		coupledIndex[nTotal]=iTempIndex;
		nExc++;
		coupledClass[nTotal]=1;
		nTotal++;
	}
	else{
		if(m_gProperty[iTempIndex]&NF_INHIBITION){
			coupledIndex[nTotal]=iTempIndex;
			nInh++;
			coupledClass[nTotal]=2;
			nTotal++;
		}
	}
}*/

void ChemicalSNetwork::DivideExcAndInh(){
	int iDimension=sqrt(m_nNeuron)+0.5;
	vector<int> exc_vec,inh_vec,nocoupling_vec;
	int TempIndex;
	for(int i=0;i<m_nNeuron;++i){
		if(m_gProperty[i]&NF_ISCOUPLED){	
			for(int j=0;j<m_gCoupleNum[i];++j){
				TempIndex=m_gCoupleIdx[i][j];
				if(m_gProperty[TempIndex]&NF_EXCITATION){
					exc_vec.push_back(TempIndex);
				}
				else if(m_gProperty[TempIndex]&NF_INHIBITION){
					inh_vec.push_back(TempIndex);
				}
				else{
					nocoupling_vec.push_back(TempIndex);
				}
			}
			aExcitation[i]=exc_vec.size();
			aInhibition[i]=inh_vec.size();
			int idx=0;
			for(int iter=0;iter<exc_vec.size();++iter){
				m_gCoupleIdx[i][idx]=exc_vec[iter];
				idx++;
			}
			for(int iter=0;iter<inh_vec.size();++iter){
				m_gCoupleIdx[i][idx]=inh_vec[iter];
				idx++;
			}
			for(int iter=0;iter<nocoupling_vec.size();++iter){
				m_gCoupleIdx[i][idx]=nocoupling_vec[iter];
				idx++;
			}
			exc_vec.clear();
			inh_vec.clear();
			nocoupling_vec.clear();
		}
		else{
			aExcitation[i]=0;
			aInhibition[i]=0;
		}
	}
}

void ChemicalSNetwork::BuildSquareSparser(const double _sparse){
	SparseNetwork::BuildSquareSparser(_sparse);
	DivideExcAndInh();
}
void ChemicalSNetwork::BuildSquareSmallWorld(const double _rewiring){
	SparseNetwork::BuildSquareSmallWorld(_rewiring);
	DivideExcAndInh();
}

void ChemicalSNetwork::BuildSquare(){
	SparseNetwork::BuildSquare();
	DivideExcAndInh();
}


//the summation being performed over all the spikes emitted by presynatic neurons at times t_spike.
//the synaptic interaction is usually classified according to whether Vsyn is larger or small than the m_threshold potential Vth,at which the postsynaptic neuron generates spikes.
//for Vsyn>Vth the interaction is called excitatory ,while for Vsyn<Vth it is called inhibitory
SpikeSNetwork::SpikeSNetwork():ChemicalSNetwork(){
	m_gSpikeTime=new double[m_nNeuron];
	memset(m_gSpikeTime,0,sizeof(double)*m_nNeuron);
	m_gVpre=new double[m_nNeuron];
	memset(m_gVpre,0,sizeof(double)*m_nNeuron);
	m_gVpre1=new double[m_nNeuron];
	memset(m_gVpre1,0,sizeof(double)*m_nNeuron);
}

SpikeSNetwork::SpikeSNetwork(int _Connor,int _HH,int _ML1,int _ML2,int _Inh,int _InhNoCoupling):ChemicalSNetwork(_Connor,_HH,_ML1,_ML2,_Inh,_InhNoCoupling){
	m_gSpikeTime=new double[m_nNeuron];
	memset(m_gSpikeTime,0,sizeof(double)*m_nNeuron);
	m_gVpre=new double[m_nNeuron];
	memset(m_gVpre,0,sizeof(double)*m_nNeuron);
	m_gVpre1=new double[m_nNeuron];
	memset(m_gVpre1,0,sizeof(double)*m_nNeuron);
}

SpikeSNetwork::~SpikeSNetwork(){
	delete []m_gSpikeTime;
	delete []m_gVpre;
	delete []m_gVpre1;
}

void SpikeSNetwork::__DoSyncVar(){
	m_t=0.0;
	for(int i=0;i<m_nNeuron;++i)
		m_gNeuron[i]->InitVar();
	srand(unsigned int(time(NULL)));
	for(int i=0;i<m_nNeuron;++i){
		m_gVpre[i]=m_gVpre1[i]=m_gNeuron[i]->V;
    }
	for(int i=0;i<m_nNeuron;++i){
		m_gSpikeTime[i]=-DBL_MAX;
	}
}

void SpikeSNetwork::__DoAlSyncVar(){
	m_t=0.0;
	srand((unsigned int)100);               //random generator seed is 100
	for(int i=0;i<m_nNeuron;++i)
		m_gNeuron[i]->RandVar();
	for(int i=0;i<m_nNeuron;++i){
		m_gVpre[i]=m_gVpre1[i]=m_gNeuron[i]->V;
    }
	for(int i=0;i<m_nNeuron;++i){
		m_gSpikeTime[i]=-DBL_MAX;
	}
//	printf("variation is initialized normalA is %lf\n",NormalA);
}

void SpikeSNetwork::UpdateSpike(){
	for(int i=0;i<m_nNeuron;++i){
		if(m_gVpre[i]>0.5&&m_gVpre[i]>m_gNeuron[i]->V&&m_gVpre[i]>m_gVpre1[i])
			m_gSpikeTime[i]=m_t;
		m_gVpre1[i]=m_gVpre[i];
		m_gVpre[i]=m_gNeuron[i]->V;
	}
}

void SpikeSNetwork::SetSpike(double _tau1,double _tau2,double _gc_exc,double _gc_inh,double _Vsyn_exc,double _Vsyn_inh){
	m_tau1=_tau1;
	m_tau2=_tau2;
	m_GcExc=_gc_exc;
	m_GcInh=_gc_inh;
	m_NormalA=1.0/(exp(-(m_tau2/(m_tau1-m_tau2))*log(m_tau1/m_tau2))-exp(-(m_tau1/(m_tau1-m_tau2))*log(m_tau1/m_tau2)));
	m_VsynExc=_Vsyn_exc;
	m_VsynInh=_Vsyn_inh;
}


void SpikeSNetwork::__DoUpdateCouple(){
	memset(sm_gCouple,0,sizeof(double)*m_nNeuron);
	UpdateSpike();
	double fCoupleExc,fCoupleInh;
	int nCount,nExcitation;
	for(int i=0;i<m_nNeuron;++i){
		if(m_gProperty[i]&NF_ISCOUPLED){
			fCoupleExc=0.0;fCoupleInh=0.0;
			nExcitation=aExcitation[i];
			nCount=aExcitation[i]+aInhibition[i];
			for(int j=0;j<nExcitation;++j){
				fCoupleExc+=m_NormalA*(exp(-(m_t-m_gSpikeTime[j])/m_tau1)-exp(-(m_t-m_gSpikeTime[j])/m_tau2));
			}
			for(int j=nExcitation;j<nCount;++j){
					fCoupleInh+=m_NormalA*(exp(-(m_t-m_gSpikeTime[j])/m_tau1)-exp(-(m_t-m_gSpikeTime[j])/m_tau2));
			}
			sm_gCouple[i]=fCoupleExc*m_GcExc*(m_VsynExc-m_gNeuron[i]->V)+fCoupleInh*m_GcInh*(m_VsynInh-m_gNeuron[i]->V);
//		printf("CoupleVector[%d] is %lf\n",i,CoupleVector[i]);
		}
	}
}

void SpikeSNetwork::__DoMakeDirectory(char *direct){
	sprintf_s(direct,100,"F:\\output\\Spike\\%s\\Connor%dHH%dML1_%dML2_%dInh_(%d,%d)",m_sName,m_nConnor,	m_nHH,m_nML1,m_nML2,m_nML2_Inh,m_nML2_Inh_Isolated);
	_mkdir("F:\\output\\Spike");
	char direct1[50];
	sprintf_s(direct1,50,"F:\\output\\Spike\\%s",m_sName);
	_mkdir(direct1);
	_mkdir(direct);
}

void SpikeSNetwork::__DoMakeFilename(char *filename){
	sprintf_s(filename,50,"gc_exc=%.5lf_gc_inh=%.5lf",m_GcExc,m_GcInh);
}

///*Update spike*/
//void Network::UpdateSpikeEx(){
//	for(int i=0;i<m_nNeuron;++i){
//		if(m_gVpre[i]>0.5&&m_gVpre[i]>m_gNeuron[i]->V&&m_gVpre[i]>m_gVpre1[i]){
//			if(spike_num[i]!=m_nNeuron){
//				m_gSpikeTime[i][spike_num[i]]=m_t;
//				spike_num[i]=spike_num[i]+1;
//			}
//			else{
//				memmove(m_gSpikeTime[i],m_gSpikeTime[i]+1,sizeof(double)*(m_nNeuron-1));
//				m_gSpikeTime[i][m_nNeuron-1]=m_t;
//			}
//		}
//		m_gVpre1[i]=m_gVpre[i];
//		m_gVpre[i]=m_gNeuron[i]->V;
//	}
//}
//
//
///*Update network synaptic coupling terms*/
//void Network::UpdateChemicalEx(){
//	memset(sm_gCouple,0,sizeof(double)*m_nNeuron);
//	int nCount;
//	double fCouple;
//	for(int i=0;i<m_nNeuron;++i){
//		nCount=aExcitation[i]+aInhibition[i];
//		fCouple=0.0;
//		for(int j=0;j<nCount;++j){
//				for(int k=0;k<spike_num[j];++k)
//					fCouple+=NormalA*(exp(-(m_t-m_gSpikeTime[j][k])/tau1)-exp(-(m_t-m_gSpikeTime[j][k])/tau2));
//		}
//		sm_gCouple[i]=fCouple*gc_c*(m_VsynExc-m_gNeuron[i]->V);
////   	printf("NetworkCouple[%d]=%lf\n",i,NetworkCouple[i]);
//	}
//}



SigmoidalSNetwork::SigmoidalSNetwork():ChemicalSNetwork(){
}
SigmoidalSNetwork::SigmoidalSNetwork(int _Connor,int _HH,int _ML1,int _ML2,int _Inh,int _InhNoCoupling):ChemicalSNetwork(_Connor,_HH,_ML1,_ML2,_Inh,_InhNoCoupling){
}

SigmoidalSNetwork::~SigmoidalSNetwork()
{
}

void SigmoidalSNetwork::SetSigmoidal(double _threshold,double _gc_exc,double _gc_inh,double _Vsyn_exc,double _Vsyn_inh){
	m_threshold=_threshold;
	m_GcExc=_gc_exc;
	m_GcInh=_gc_inh;
	m_VsynExc=_Vsyn_exc;
	m_VsynInh=_Vsyn_inh;
}

void SigmoidalSNetwork::__DoUpdateCouple(){
	memset(sm_gCouple,0,sizeof(double)*m_nNeuron);
	double fCoupleExc,fCoupleInh;
	int nExc,nTotal,iIndex;
	for(int i=0;i<m_nNeuron;++i){
		if(m_gProperty[i]&NF_ISCOUPLED){
			fCoupleExc=0.0;fCoupleInh=0.0;
			nExc=aExcitation[i];nTotal=aExcitation[i]+aInhibition[i];
			for(int j=0;j<nExc;++j){
				iIndex=m_gCoupleIdx[i][j];
				fCoupleExc+=1.0/(1.0+exp(-(m_gNeuron[iIndex]->V-m_threshold)));
			}
			for(int j=nExc;j<nTotal;++j){
				iIndex=m_gCoupleIdx[i][j];
				fCoupleInh+=1.0/(1.0+exp(-(m_gNeuron[iIndex]->V-m_threshold)));
			}
			sm_gCouple[i]=fCoupleExc*m_GcExc*(m_VsynExc-m_gNeuron[i]->V)+fCoupleInh*m_GcInh*(m_VsynInh-m_gNeuron[i]->V);
		}
	}
}

void SigmoidalSNetwork::__DoMakeDirectory(char *direct){
	sprintf_s(direct,100,"F:\\output\\Sigmoidal\\%s\\Connor%dHH%dML1_%dML2_%dInh_(%d,%d)",m_sName,m_nConnor,m_nHH,m_nML1,m_nML2,m_nML2_Inh,m_nML2_Inh_Isolated);
	_mkdir("F:\\output\\Sigmoidal");
	char direct1[50];
	sprintf_s(direct1,"F:\\output\\Sigmoidal\\%s",m_sName);
	_mkdir(direct1);
	_mkdir(direct);
}

void SigmoidalSNetwork::__DoMakeFilename(char *filename){
	sprintf_s(filename,40,"gc_exc=%.5lf_gc_inh=%.5lf",m_GcExc,m_GcInh);
}
