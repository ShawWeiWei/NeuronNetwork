#pragma once
#include "Network.h"

class ChemicalSNetwork:public SparseNetwork{
public:
	ChemicalSNetwork();
	ChemicalSNetwork(int _Connor,int _HH,int _ML1,int _ML2,int _Inh,int _InhNoCoupling);
	virtual ~ChemicalSNetwork();
	void BuildSquare();
	void BuildSquareSmallWorld(const double _rewiring);
	void BuildSquareSparser(const double _sparse);
//	void ModifyConnection();
protected:
	int *aExcitation;
	int *aInhibition;	
private:
	void DivideExcAndInh();
};

class SpikeSNetwork:public ChemicalSNetwork{
public:
	SpikeSNetwork();
	SpikeSNetwork(int _Connor,int _HH,int _ML1,int _ML2,int _Inh,int _InhNoCoupling);
	virtual ~SpikeSNetwork();
	void SetSpike(double _tau1,double _tau2,double _gc_exc,double _gc_inh,double _Vsyn_exc,double _Vsyn_inh);
private:
	//the summation being performed over all the spikes emitted by the presynaptic neurons
	double m_tau1;
	double m_tau2;
	double m_NormalA;
	double m_GcExc;
	double m_GcInh;
	double m_VsynExc;
	double m_VsynInh;

//	int*spike_num;
	double *m_gSpikeTime;
	double *m_gVpre;
	double *m_gVpre1;

	virtual void __DoUpdateCouple();
	void UpdateSpike();
	virtual void __DoSyncVar();
	virtual void __DoAlSyncVar();
	virtual void __DoMakeDirectory(char *direct);
	virtual void __DoMakeFilename(char *filename);
};


class SigmoidalSNetwork:public ChemicalSNetwork{
public:
	SigmoidalSNetwork();
	SigmoidalSNetwork(int _Connor,int _HH,int _ML1,int _ML2,int _Inh,int _InhNoCoupling);
	virtual ~SigmoidalSNetwork();
	void SetSigmoidal(double _threshold,double _gc_exc,double _gc_inh,double _Vsyn_exc,double _Vsyn_inh);
private:
	//sigmoidal
	double m_threshold;
	double m_GcExc;
	double m_GcInh;
	double m_VsynExc;
	double m_VsynInh;

private:
	virtual void __DoUpdateCouple();
	virtual void __DoMakeDirectory(char *direct);
	virtual void __DoMakeFilename(char *filename);
};