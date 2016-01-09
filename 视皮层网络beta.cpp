// 视皮层网络.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "Network.h"
#include "Neuron.h"
#include "MorrisLecarModel.h"
#include "HHModel.h"
#include "ConnorModel.h"
#include "ElectricalNetwork.h"
#include "ChemicalNetwork.h"


/*test function*/
void TestMorrisLecar1();
void TestMorrisLecar2();
void TestConnor();
void TestHH();
void TestSquare();
void TestNetwork();
void TestRandShuffle();
void TestNoise();
void TestEfficiency();
/*Square Heterogenous*/
void SquareHeterogenous();
void DynamicalNoiseIntensity(const double _gc_s);


int _tmain(int argc, _TCHAR* argv[])
{
	double t_begin,t_end;
	t_begin=clock();	

/*	ElectricalSNetwork instance(0,0,16384,0);
	instance.BuildSquare();
	instance.SetElectrical(0.85);
	instance.SetNoiseIntensity(1.5);
	instance.SpiralWave();*/
	

	
	SigmoidalSNetwork instance(0,0,15565,0,819,0);// 16384--->100% 15565--->95% 14746--->90% 13107--->80% 3277--->20% 1638 ---->10% 819---->5% 655--->4% 492--->3% 328--->2% 164--->1%
//	instance.BuildSquare();
	instance.BuildSquareSmallWorld();
	instance.OutputCoupleIdx();
//	instance.SetNoiseIntensity(0.00);
//	instance.SetSigmoidal(-25,0.4,0.4,30,-60);
//	instance.SpiralWave();
	//for(int i=0;i<10;++i){
	//	if(4==i) continue;
	//	instance.SetSigmoidal(-25,0.16+i*0.01,0.16+i*0.01,30,-60);
	//	instance.SpiralWave();
	//}

//	instance.MaximumPotential();
//    TestNetwork();
//	TestRandShuffle();
//	SquareHeterogenous();
//	TestSquare();
//	Network homo_network(0,0,15565,0,164,655);// 16384--->100% 15565--->95% 819---->5% 655--->4% 492--->3% 328--->2% 164--->1%
//	homo_network.BuildSquare();
	
//	homo_network.OutputSquare();
//	homo_network.SetElectrical(0.75);
//	homo_network.SetChemical(2.0,1.0,1.0);
//	homo_network.SetSigmoidal(-25,0.265,30,-60);
//    homo_network.SetNoiseIntensity(0.002);
//	homo_network.SetDt(0.01);

//    homo_network.SpiralWaveAndTimeSeries();
//	homo_network.OutputTimeSeries();
//	MorrisLecar mlneuron;
//	mlneuron.set_class1();
//	mlneuron.OutputTimeSeries(39.97);
//	HH neuronmodel;
	//neuronmodel.set_class1();
//	neuronmodel.Frequency_I(8.0,12.0,0.01);
//	TestNoise();
//	double fCoupleIntensity[]={0.24,0.25,0.255,0.26,0.265,0.27,0.275};
//	int nCount=sizeof(fCoupleIntensity)/sizeof(double);
//	printf("nCount is %d\n",nCount);
//	for(int i=0;i<nCount;++i){
//		DynamicalNoiseIntensity(fCoupleIntensity[i]);
//		printf("%d is done!\n",i);
//	}
	t_end=clock();
	printf("Process exited after %.3lf seconds\n",(t_end-t_begin)/1000);                      
	system("PAUSE");
	return 0;
}

void TestMorrisLecar1(){
	MorrisLecar ML1;
	ML1.set_class1();
	ML1.Frequency_I(35,45,0.01);
}
void TestMorrisLecar2(){
	MorrisLecar ML2;
	ML2.set_class2();
	ML2.Frequency_I(85,90,0.01);

}
void TestConnor(){
	Connor con;
	con.Frequency_I(0,40,0.1);
}
void TestHH(){
	HH hh;
	hh.Frequency_I(5,15,0.01);
	hh.OutputTimeSeries(0);
}

//void TestSquare(){
//	Network Square_network(0,0,50,50,0,0);
//	Square_network.BuildSquare();
//	Square_network.OutputNetworkCouple();
//}
//
//void TestNetwork(){
//	Network Square_ML1_network(0,0,100,0,0,0);
//	Square_ML1_network.OutputTimeSeries();
//}

void TestRandShuffle(){
	int array1[20];
	for(int i=0;i<20;++i)
		array1[i]=i;
	std::random_shuffle(array1,array1+20);
	for(int i=0;i<20;++i)
		printf("%d ",array1[i]);
}

void TestNoise(){
	int number=100000;
	double *noise=new double[number];
	memset(noise,0,sizeof(double)*number);
	double ran1,ran2;
	srand((unsigned int)1000);
	for(int i=0;i<number/2;++i){
		ran1=Network::Uniform_01();
		ran2=Network::Uniform_01();
		noise[i]=10.0+sqrt(10.0)*sqrt(-2*log(ran1))*cos(2*PI*ran2);
		noise[i+number/2]=10.0+sqrt(10.0)*sqrt(-2*log(ran1))*sin(2*PI*ran2);
	}
	if(number&0x1){
		ran1=Network::Uniform_01();
		ran2=Network::Uniform_01();
		noise[number-1]=10.0+sqrt(10.0)*sqrt(-2*log(ran1))*cos(2*PI*ran2);
	}

	FILE *noisefile;
	fopen_s(&noisefile,"noise_file.dat","w");
	for(int i=0;i<number;++i)
		fprintf(noisefile,"%lf\n",noise[i]);
	fclose(noisefile);

	delete []noise;
}

//void testefficiency(){
//	ElectricalSNetwork homo_network(0,0,15565,0);// 16384--->100% 15565--->95% 819---->5% 655--->4% 492--->3% 328--->2% 164--->1%
//	homo_network.BuildSquare();
//	for(int i=0;i<100000;i++)
//		homo_network.UpdateNoise();
//}
void SquareHeterogenous(){
	SpikeSNetwork SH_network(0,0,5000,5000,0,0);
	SH_network.BuildSquare();
	SH_network.SetSpike(3.0,1.0,0.01,0.2,30,-60);
	SH_network.SpiralWave();
	SH_network.SetSpike(3.0,1.0,0.05,0.2,30,-60);
	SH_network.SpiralWave();
}

void DynamicalNoiseIntensity(const double _gc_s){
	SigmoidalSNetwork homo_network(0,0,16384,0,0,0);
	homo_network.BuildSquare();
	homo_network.SetSigmoidal(-25,_gc_s,10.0,30,-60);
	double order;
	for(int i=0;i<5;++i){
		order=pow(10,-i);
		homo_network.SetNoiseIntensity(order);
		homo_network.SpiralWaveAndTimeSeries();
		homo_network.SetNoiseIntensity(2.0*order);
		homo_network.SpiralWaveAndTimeSeries();
		homo_network.SetNoiseIntensity(5.0*order);
		homo_network.SpiralWaveAndTimeSeries();
	}
}