#pragma once
#include <stdlib.h>
#include <memory.h>
#include <float.h>
#include <direct.h>
#include <math.h>
#include <time.h>
#include "Neuron.h"
#include <vector>
#include <algorithm>
#include <string.h>
using std::vector;
typedef unsigned char u8;
#define NF_CLASS_ONE 0x01
#define NF_CLASS_TWO 0x02
#define NF_RESTING 0x04
#define NF_FIRING 0x08
#define NF_ISCOUPLED 0x10
#define NF_EXCITATION 0x20
#define NF_INHIBITION 0x40
#define NF_NOCOUPLING 0x80


struct networkProperty{
	unsigned networkAttribute;
	int nMl1;
	int nMl2;
	int nInh;
	int nInhNoCoupling;
};

#define PI 3.1415926535898

class neuron;
class Network{
	friend class Connor;
	friend class HH;
	friend class MorrisLecar;
public:
	/*Default Constructor*/
	Network();
	Network(const int _Connor,const int _HH,const int _ML1,const int _ML2,const int _Inh,const int _InhNoCoupling);
	/*Destructor*/
	virtual ~Network();       
	//Preprocess 2(optional):Specify synaptical decay time, rise time and coupling intensity
	void SetNoiseIntensity(const double _noist_intensity);
	//Preprocess 3:Specify neural currents
	void ConfigureI(double _ML1_I_Span,double _ML2_I_Span);
	void RandI();
	/*Set time step length...*/
	void SetDt(const double&_dt=0.01);
	//Output network for verify
	void OutputCouple();
	void OutputNoise();


	void OutputNoForOneAndTwo();
	static double Uniform_01(){
		return ((rand()+1.0)/(RAND_MAX+1.0));
	}

	//process
	void OutputI();
	void OutputTimeSeries();
	void SpiralWave(int _begin,int _end,int _dt);
	void SpiralWave();
	void SpiralWaveAndTimeSeries();
	void OutputSpikingIndex();
	void OutputPhaseAmplitude();
	void OutputAverISI();
	void OutputPopulationFirings();
	void OutputPopulationFiringsOnce();
	void MaximumPotential();
	void Tau1Tau2SyncFactor(const double _tau_max,const double _interval,const double _g);
	double SyncFactor();
	double OrderParameter(const double _eval_begin,const double _iter_end);
	void Tau1Tau2OrderParameter(const double _tau_max,const double _interval,const double _g);
	void OutputISI();

protected:
	const int m_nConnor;
	const int m_nHH;
	const int m_nML1;
	const int m_nML2;
	const int m_nML2_Inh;
	const int m_nML2_Inh_Isolated;
	const int m_nNeuron;

	double m_t;
	static double sm_dt;

	double m_noiseintensity;
	double m_coupleintensity;


	neuron **m_gNeuron;
	static double *sm_gCouple;
	static double *sm_gNoise;
	char *m_sName;
	u8 *m_gProperty;

		//Config neurons' currents
	bool m_isRandI;
	double m_ML1_I_Begin;
	double m_ML1_I_Span;
	double m_ML2_I_Begin;
	double m_ML2_I_Span;	
private:

	//Network  iteration
	void __EulerIterate();
	void __RungekuttaIterate();
	void __UpdateNoise();

	bool __IsML1(int No);
	bool __IsML2(int No);
	bool __IsInner(int No);
	//Preprocess 4:Set the initial values of variables to sychronization...
	void __DoSyncVar();
	void __DoAlSyncVar();
	virtual void __DoUpdateCouple()=0;
	virtual void __DoMakeDirectory(char *direct)=0;
	virtual void __DoMakeFilename(char *filename)=0;
};


class SparseNetwork:public Network{
protected:
	int ** m_gCoupleIdx;
	int *m_gCoupleNum;

public:
	SparseNetwork();
	SparseNetwork(int _HH,int _Connor,int _ML1,int _ML2,int _Inh,int _InhNoCoupling);
	virtual ~SparseNetwork();
	void BuildNearestNeighbor(const int _adjacent);  
	void BuildSquareSmallWorld(const double _rewiring);   
	void BuildSquareSparser(const double _double);
	void BuildSquare();
	void OutputCoupleIdx();


private:
//	virtual void __DoUpdateCouple()=0;
};

