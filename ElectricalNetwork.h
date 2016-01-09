#pragma once

#include "Network.h"

class ElectricalSNetwork:public SparseNetwork{
public:
	ElectricalSNetwork();
	ElectricalSNetwork(int _Connor,int _HH,int _ML1,int _ML2);
	void SetElectrical(const double _gc_e);

private:
		//electrical
//	double gc_e;

private:
	virtual void UpdateCouple();
	virtual void contextualDirectory(char *direct);
	virtual void SpecifyFilename(char *filename);
};