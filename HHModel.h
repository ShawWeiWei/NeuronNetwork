#pragma once

#include "Neuron.h"

/*invariant parameters*/
/*#define V_NA 50.0
#define V_K -77.0
#define V_L -54.4
#define G_NA 120.0
#define G_K 36.0
#define G_L 0.3
#define C 1.0*/

class HH:public neuron{
	friend class Network;
public:
	HH();
	HH(const int&_No,bool _IsExcitation);
	~HH();
	void InitVar();
	void RandVar();
	void equation_s();
	void equation();
	void euler_s();
	void euler();
	void rungekutta_s();
	void rungekutta();
	double alpha_m();
	double beta_m();
	double alpha_h();
	double beta_h();
	double alpha_n();
	double beta_n();
private:
	const int No;


	double m;
	double h;
	double n;
	double DV;
	double Dm;
	double Dh;
	double Dn;
	double V0;
	double m0;
	double h0;
	double n0;
	double DV1;
	double Dm1;
	double Dh1;
	double Dn1;
	double DV2;
	double Dm2;
	double Dh2;
	double Dn2;
	double DV3;
	double Dm3;
	double Dh3;
	double Dn3;
	double DV4;
	double Dm4;
	double Dh4;
	double Dn4;
};

inline double HH::alpha_m(){
	return 0.1*(V+40.0)/(1.0-exp((-V-40.0)/10.0));
}

inline double HH::beta_m(){
	return 4.0*exp((-V-65.0)/18.0);
}

inline double HH::alpha_h(){
	return 0.07*exp((-V-65.0)/20.0);
}

inline double HH::beta_h(){
	return 1.0/(1.0+exp((-V-35.0)/10.0));
}

inline double HH::alpha_n(){
	return 0.01*(V+55.0)/(1.0-exp((-V-55.0)/10.0));
}

inline double HH::beta_n(){
	return 0.125*exp((-V-65.0)/80.0);
}