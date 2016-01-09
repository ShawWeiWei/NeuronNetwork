#pragma once
#include "Neuron.h"


/*Invariant PARAMETERS in the ML model*/
/*#define C 20.0//C is the membrane capacitance
#define gK 8.0//gK is maximum potassium conductance
#define gL 2.0//gL is leak conductance
#define VK -84.0//VK is reversal potential for the potassium current
#define VCa 120.0//VCa is reversal potential for calcium current
#define VL -60.0//VL is membrane resting potential 
#define V1 -1.2
#define V2 18.0*/



class MorrisLecar:public neuron{
	friend class Network;
public:
	MorrisLecar(); 
	MorrisLecar(const int&_No,bool _IsClass1,bool _IsResting);
	~MorrisLecar();
	void equation();
	void euler();
	void rungekutta();

//	void equation_s();
	void euler_s();
	void rungekutta_s();

//	void set_dt(const double&_dt=0.01);
    void InitVar();
	void RandVar();
//	void ResetNoise();
	void set_class1();
	void set_class2();

private:

	const int No;


	double n;
	/*parameters*/
	double gCa;
	double phi;
	double V3;
	double V4;


	/*intermediate variables for equations*/
	double DV;
	double Dn;
	
	/*intermediate variables for updating V and n*/
	double V0;
	double n0;

	/*intermediate variables for rungekutta method*/
	double DV1;
	double Dn1;
	double DV2;
	double Dn2;
	double DV3;
	double Dn3;
	double DV4;
	double Dn4;

	/*variables to compare to determine spike time*/
	/*double V_pre;
	double V_pre1;*/
};