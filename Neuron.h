#pragma once
#include "math.h"
#include "string.h"
#include "Network.h"

class neuron{
public:
	neuron();
	virtual ~neuron();

	virtual void euler_s()=0;
	virtual void euler()=0;
	virtual void rungekutta_s()=0;
	virtual void rungekutta()=0;
	virtual void InitVar()=0;
	virtual void RandVar()=0;

	void SetI(double _I);
	void SetRandI(double _rand,double _begin,double _span);
	void OutputTimeSeries(const double&_I);
	void Frequency_I(const double&_I_Begin,const double&_I_End,const double&_I_Partition);
	void OutputISI();
	void OutputISIs(const double&_I_Begin,const double&_I_End,const double&_I_Partition);
//private:
//	u8 eProperty;
	double V;
	double I;
	char neuron_model[15];
protected:
	
	const double single_dt;
	

};

