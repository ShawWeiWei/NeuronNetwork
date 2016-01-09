#pragma once

#include "Neuron.h"


/*invariant parameters*/
//#define V_NA 55.0
//#define V_K -72.0
//#define V_L -17.0
//#define V_A -75.0
//#define G_NA /*1.2*/ 120.0
//#define G_K /*0.2*/ 20.0
//#define G_L /*0.003*/ 0.3
//#define G_A /*0.477*/ 47.7
//#define C /*10.0*/ 1.0




class Connor:public neuron{
	friend class Network;
public:
	Connor();
	Connor(const int&_No,bool _IsExcitation);
	~Connor();
	virtual void InitVar();
	virtual void RandVar();
	virtual void equation_s();
	virtual void equation();
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
	double a_infinity();
	double tau_a();
	double b_infinity();
	double tau_b();
private:
	const int No;

	double m;
	double h;
	double n;
	double A;
	double B;
	double DV;
	double Dm;
	double Dh;
	double Dn;
	double DA;
	double DB;
	double V0;
	double m0;
	double h0;
	double n0;
	double A0;
	double B0;
	double DV1;
	double Dm1;
	double Dh1;
	double Dn1;
	double DA1;
	double DB1;
	double DV2;
	double Dm2;
	double Dh2;
	double Dn2;
	double DA2;
	double DB2;
	double DV3;
	double Dm3;
	double Dh3;
	double Dn3;
	double DA3;
	double DB3;
	double DV4;
	double Dm4;
	double Dh4;
	double Dn4;
	double DA4;
	double DB4;
};

inline double Connor::alpha_m(){
	return 0.1*(V+29.7)/(1.0-exp((-V-29.7)/10.0));
}

inline double Connor::beta_m(){
	return 4.0*exp((-V-54.7)/18.0);
}

inline double Connor::alpha_h(){
	return 0.07*exp((-V-48.0)/20.0);
}

inline double Connor::beta_h(){
	return 1.0/(1.0+exp((-V-18.0)/10.0));
}

inline double Connor::alpha_n(){
	return 0.01*(V+46.7)/(1.0-exp((-V-46.7)/10.0));
}

inline double Connor::beta_n(){
	return 0.125*exp((-V-56.7)/80.0);
}

inline double Connor::a_infinity(){
	return pow(0.0761*exp((V+94.22)/31.84)/(1.0+exp((V+1.17)/28.93)),1.0/3.0);
}

inline double Connor::tau_a(){
	return 0.3632+1.158/(1.0+exp((V+55.96)/20.12));
}

inline double Connor::b_infinity(){
	return 1.0/pow(1.0+exp((V+53.3)/14.54),4.0);
}

inline double Connor::tau_b(){
	return 1.24+2.678/(1.0+exp((V+50.0)/16.027));
}