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
void Configure(SigmoidalSNetwork *pInstance,double coupleintensity,double noise,double i_span){
	pInstance->SetNoiseIntensity(noise);
	pInstance->SetSigmoidal(-25,coupleintensity,coupleintensity,30,-60);
	pInstance->BuildSquare();
	pInstance->ConfigureI(i_span,i_span);
	pInstance->RandI();
}
void HeterNoNoise();
void I_NoNoise();
int _tmain(int argc, _TCHAR* argv[])
{
	double t_begin,t_end;
	t_begin=clock();	

	SigmoidalSNetwork instance(0,0,16384,0,0,0);
	instance.BuildSquare();
	instance.SetNoiseIntensity(0.0);
	instance.ConfigureI(0.0,0.0);
	instance.RandI();
	instance.SetSigmoidal(-25,0.4,0.4,30,-60);
	instance.OutputCouple();

	t_end=clock();
	printf("Process exited after %.3lf seconds\n",(t_end-t_begin)/1000);                      
	system("PAUSE");
	return 0;
}

void HeterNoNoise(){
	double aGc[]={0.1,0.2,0.3,0.4,0.5};//{0.1,0.2,0.3,0.4,0.5};
	double aI[]={0.2};
	int aML2[]={1,2,3,4};
	int sizeGc=sizeof(aGc)/sizeof(double);
	int sizeI=sizeof(aI)/sizeof(double);
	int sizeML2=sizeof(aML2)/sizeof(int);

	for(int iML2=0;iML2<sizeML2;iML2++){
		int nML2=16384.0*aML2[iML2]/100.0+0.5;
		int nML1=16384-nML2;
		SigmoidalSNetwork instance(0,0,nML1,nML2,0,0);
	//	instance.OutputNoForOneAndTwo();
		instance.BuildSquare();
		instance.SetNoiseIntensity(0.0);
		for(int iI=0;iI<sizeI;++iI){
			instance.ConfigureI(aI[iI],aI[iI]);
			instance.RandI();
	//		instance.OutputI();
			for(int iGc=0;iGc<sizeGc;++iGc){
				instance.SetSigmoidal(-25,aGc[iGc],aGc[iGc],30,-60);
				instance.OutputNoForOneAndTwo();
				instance.OutputTimeSeries();
//				instance.SpiralWaveAndTimeSeries();
//				instance.OutputSpikingIndex();
//				instance.OutputPhaseAmplitude();
//				instance.OutputAverISI();
	//			instance.OutputPopulationFiringsOnce();
	//			instance.SpiralWaveAndTimeSeries();
			}
		}
	}
}

void I_NoNoise(){
	double aGc[]={0.1,0.2,0.3,0.4,0.5};
	int sizeOfaGc=sizeof(aGc)/sizeof(double);
	double aI[]={0.0,0.1,0.2,0.3};//{0.2,0.3};
	int sizeOfaI=sizeof(aI)/sizeof(double);
	for(int prop=2;prop<5;++prop){
		int nML2=16384.0*prop/100.0+0.5;
		int nML1=16384-nML2;
		SigmoidalSNetwork instance(0,0,nML1,nML2,0,0);
		instance.BuildSquare();
		instance.SetNoiseIntensity(0.0);
		for(int i=0;i<sizeOfaI;++i){
			instance.ConfigureI(aI[i],aI[i]);
			instance.RandI();
			for(int j=0;j<sizeOfaGc;++j){
			instance.SetSigmoidal(-25,aGc[j],aGc[j],30,-60);
			instance.SpiralWaveAndTimeSeries();
			instance.OutputSpikingIndex();
			instance.OutputPhaseAmplitude();
			instance.OutputAverISI();
			}
		}
    }
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