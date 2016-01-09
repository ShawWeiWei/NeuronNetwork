#include "stdafx.h"
#include "Neuron.h"

neuron::neuron():single_dt(0.01){
}

neuron::~neuron(){}

void neuron::SetI(double _I){
	I=_I;
}

void neuron::SetRandI(double _rand,double _begin,double _span){
	I=_begin+_span*_rand;
//	printf("%lf\n",I);
}

void neuron::OutputTimeSeries(const double&_I){
    /*Initialize...*/
	InitVar();
	I=_I;


	/*Create directory...*/
	/*char dir[30];
	sprintf_s(dir,"I=%.5lfmV",_I);*/
	_mkdir(neuron_model);

	char Outfilename[50];
	sprintf_s(Outfilename,".\\%s\\I=%.5lfmV_TimeSeries.dat",neuron_model,_I);
	FILE *pOutputTimeSeries;
	fopen_s(&pOutputTimeSeries,Outfilename,"w");
	fprintf(pOutputTimeSeries,"%lf %lf\n",0,V);
	for(int i=1;i<600000;++i){
		rungekutta_s();
		fprintf(pOutputTimeSeries,"%lf %lf\n",i*single_dt,V);
	}
	fclose(pOutputTimeSeries);
}

void neuron::Frequency_I(const double&_I_Begin,const double&_I_End,const double&_I_Partition){
	/*Initialize files to save ISIs.......*/
	FILE *ISI;
	_mkdir(neuron_model);
	char Outfilename[50];
	sprintf_s(Outfilename,".\\%s\\Frequency-I.dat",neuron_model);
	fopen_s(&ISI,Outfilename,"w");
	double V_pre,V_pre1;
	double t_f,freq;
	double max,min;

	

	const int step=(int)((_I_End-_I_Begin)/_I_Partition);
	for(int i=0;i<step;++i){
		/*Initialize.......*/
		InitVar();
		I=_I_Begin+_I_Partition*i;

		/*Skip the transient state......*/
		for(int j=0;j<100000;++j)
			euler_s();
		/*Distinguish Now!......*/
		max=min=V;
		for(int j=0;j<100000;++j){
			euler_s();
			if(max<V) max=V;
			if(min>V) min=V;
		}
		if(abs(max-min)<0.0001||max<-10.0){
			fprintf(ISI,"%lf %lf\n",I,0.0);
			printf("I=%lf f=0.0. V=%lfmV\n",I,max);
			continue;
		}
		else{
			 /*Define temporary variables to compare*/
			V_pre1=V;
			euler_s();
			V_pre=V;

			/*Define firing time*/
			t_f=0.0;
			freq=0.0;



			/*Start evaluating and saving ISIs.....*/
			for(int j=0;j<1000000;++j){
				rungekutta_s();
				if(V_pre>20.0&&V_pre>V_pre1&&V_pre>V){
					if(t_f<0.00001) t_f=j*single_dt;
					else{ 
						freq=1000/(j*single_dt-t_f);
						fprintf(ISI,"%lf %lf\n",I,freq);
						printf("I=%lf f=%lf  V=(%lf->%lf)mV.\n",I,freq,min,max);
						break;
					}
				}
				V_pre1=V_pre;
				V_pre=V;
			}
			//printf("I=%lf is evaluated.error\n",I);
		}
	}
	fclose(ISI);
}


void neuron::OutputISI(){
	FILE *ISI;
	_mkdir(neuron_model);
	char Outfilename[50];
	sprintf_s(Outfilename,".\\%s\\I=%lf_ISI.dat",I,neuron_model);
	fopen_s(&ISI,Outfilename,"w");
	fprintf(ISI,"I=%lf\n",I);
	InitVar();

	double V_pre,V_pre1;
	V_pre=V;
	rungekutta_s();
	V_pre1=V;
	for(int iter=1;iter<100000;++iter){
		rungekutta_s();
		if(V_pre>0.5&&V_pre>V&&V_pre>V_pre1){
			break;
		}
		V_pre1=V_pre;
		V_pre=V;
	}
	int time_interval=0;
	for(int iter=1;iter<1000000;++iter){
		rungekutta_s();
		++time_interval;
		if(V_pre>0.5&&V_pre>V&&V_pre>V_pre1){
			fprintf(ISI,"%lf\n",time_interval*single_dt);
			time_interval=0;
		}
		V_pre1=V_pre;
		V_pre=V;
	}
	fclose(ISI);
}

void neuron::OutputISIs(const double&_I_Begin,const double&_I_End,const double&_I_Partition){
	FILE *ISI;
	_mkdir(neuron_model);
	char Outfilename[50];
	sprintf_s(Outfilename,".\\%s\\I=%.3lf:%.3lf:%.3lf_ISIs.dat",I,neuron_model);
	fopen_s(&ISI,Outfilename,"w");

	InitVar();

	double V_pre,V_pre1;
	V_pre=V;
	rungekutta_s();
	V_pre1=V;
	for(int iter=1;iter<100000;++iter){
		rungekutta_s();
		if(V_pre>0.5&&V_pre>V&&V_pre>V_pre1){
			break;
		}
		V_pre1=V_pre;
		V_pre=V;
	}
	int time_interval=0;
	for(int iter=1;iter<1000000;++iter){
		rungekutta_s();
		++time_interval;
		if(V_pre>0.5&&V_pre>V&&V_pre>V_pre1){
			fprintf(ISI,"%lf\n",time_interval*single_dt);
			time_interval=0;
		}
		V_pre1=V_pre;
		V_pre=V;
	}
	fclose(ISI);
}