#include "stdafx.h"
#include "postproc.h"
#include "Network.h"


void OutputSpatioTemporal(const pInfoNetwork pNetwork){
	/*Create directory...*/
	char dir[50];
	sprintf_s(dir,"Connor%dHH%dML1_%dML2_%d",pNetwork->Connor_num,pNetwork->HH_num,pNetwork->ML1_num,pNetwork->ML2_num);
	_mkdir(dir);

	/*Create file pointer...*/
	FILE *pInputTimeSeries;
	FILE *pOutputSpatioTemporal;
	char Infilename[80];
	char Outfilename[80];
	sprintf_s(Infilename,".\\%s\\tau1=%.3lftau2=%.3lfg=%.3lfTimeSeries.dat",dir,pNetwork->tau1,pNetwork->tau2,pNetwork->g);
	sprintf_s(Outfilename,".\\%s\\tau1=%.3lftau2=%.3lfg=%.3lfSpatioTemporal.dat",dir,pNetwork->tau1,pNetwork->tau2,pNetwork->g);
	fopen_s(&pInputTimeSeries,Infilename,"r");
	fopen_s(&pOutputSpatioTemporal,Outfilename,"w");

	int neuron_num=pNetwork->Connor_num+pNetwork->HH_num+pNetwork->ML1_num+pNetwork->ML2_num;
	double *V=new double[neuron_num];
	double *V_pre=new double[neuron_num];
	double *V_pre1=new double[neuron_num];
	double t;
	while(fscanf_s(pInputTimeSeries,"%lf",&t)!=EOF){
		for(int i=0;i<neuron_num;++i){
			fscanf_s(pInputTimeSeries," %lf",&V[i]);
			if(V_pre[i]>0.5&&V_pre[i]>V[i]&&V_pre[i]>V_pre1[i])
				fprintf(pOutputSpatioTemporal,"%lf %d\n",t,i);
			V_pre1[i]=V_pre[i];
			V_pre[i]=V[i];
		}
        fscanf_s(pInputTimeSeries,"\n");
    }

	delete []V;
	delete []V_pre;
	delete []V_pre1;
	fclose(pInputTimeSeries);
	fclose(pOutputSpatioTemporal);
}

void OutputSpikeTime(const pInfoNetwork pNetwork){
	/*Create directory...*/
	char dir[50];
	sprintf_s(dir,"Connor%dHH%dML1_%dML2_%d",pNetwork->Connor_num,pNetwork->HH_num,pNetwork->ML1_num,pNetwork->ML2_num);
	_mkdir(dir);

	/*Create file pointer...*/
	char Infilename[80];
	sprintf_s(Infilename,".\\%s\\tau1=%.3lftau2=%.3lfg=%.3lfSpatioTemporal.dat",dir,pNetwork->tau1,pNetwork->tau2,pNetwork->g);
	FILE *pInputSpatioTemporal;
	if(0!=fopen_s(&pInputSpatioTemporal,Infilename,"r")){
		OutputSpatioTemporal(pNetwork);
		fopen_s(&pInputSpatioTemporal,Infilename,"r");
	}
	int neuron_num=pNetwork->Connor_num+pNetwork->HH_num+pNetwork->ML1_num+pNetwork->ML2_num;
	FILE **pTransientSpikeTime=new FILE*[neuron_num];
	char filename[80];
	for(int i=0;i<neuron_num;++i){
		sprintf_s(filename,".\\%s\\SpikeTime_%d.dat",dir,i);
		fopen_s(&pTransientSpikeTime[i],filename,"w+");
	}
	double time_in;
	int index_in;

	while(fscanf_s(pInputSpatioTemporal,"%lf %d\n",&time_in,&index_in)!=EOF){
		fprintf(pTransientSpikeTime[index_in],"%lf ",time_in);
    }

	fclose(pInputSpatioTemporal);
	
	FILE *pOutputSpikeTime;
	char Outfilename[80];
	sprintf_s(Outfilename,"\\%s\\tau1=%.3lftau2=%.3lfg=%.3lfSpikeTime.dat",dir,pNetwork->tau1,pNetwork->tau2,pNetwork->g);
	fopen_s(&pOutputSpikeTime,Outfilename,"w");

	for(int i=0;i<neuron_num;++i){
		while(fscanf_s(pTransientSpikeTime[i],"%lf ",&time_in)!=EOF)
			fprintf(pOutputSpikeTime,"%lf ",time_in);
		fprintf(pOutputSpikeTime,"\n");
	}

	for(int i=0;i<neuron_num;++i){
		sprintf_s(filename,".\\%s\\SpikeTime_%d.dat",dir,i);
		remove(filename);
	}
	

	delete [] pTransientSpikeTime;
}


void OutputFrequency(const pInfoNetwork pNetwork){
	//Create directory...
	char dir[50];
	sprintf_s(dir,"Connor%dHH%dML1_%dML2_%d",pNetwork->Connor_num,pNetwork->HH_num,pNetwork->ML1_num,pNetwork->ML2_num);
	_mkdir(dir);

	//Create file pointers...
	char Infilename[80];
	sprintf_s(Infilename,".\\%s\\tau1=%.3lftau2=%.3lfg=%.3lfSpikeTime.dat",dir,pNetwork->tau1,pNetwork->tau2,pNetwork->g);
	FILE *pInputSpikeTime;
	if(0!=fopen_s(&pInputSpikeTime,Infilename,"r")){
		OutputSpikeTime(pNetwork);
		fopen_s(&pInputSpikeTime,Infilename,"r");
	}

	int neuron_num=pNetwork->Connor_num+pNetwork->HH_num+pNetwork->ML1_num+pNetwork->ML2_num;
	
	double *freq=new double[neuron_num];//Dynamically allocate memory
	double t_spike;
	double t_spike_1;
	int S;//S indicate that the neuron has fired S+1 times
	char c;
	for(int i=0;i<neuron_num;++i){
		S=0;
		
		while(1){
			fpos_t pos;
			fgetpos(pInputSpikeTime,&pos);
			while((c=fgetc(pInputSpikeTime))==' ');
			if(c=='\n') break;
			else fsetpos(pInputSpikeTime,&pos);
			if(fscanf_s(pInputSpikeTime,"%lf",&t_spike)==EOF){
				break;
			}
			else{
				if(S==0) t_spike_1=t_spike;
				S++;
			}
		}
		if(S>1) freq[i]=(S-1)/(t_spike-t_spike_1);
		else freq[i]=DBL_MAX;
	}

	double networkfreq=0.0;
	for(int i=0;i<neuron_num;++i)
		networkfreq+=freq[i];
	networkfreq/=neuron_num;

	FILE *pOutputFrequency;
	char outfilename[80];
	sprintf_s(outfilename,".\\%s\\tau1=%.3lftau2=%.3lfg=%.3lfFrequency.dat",dir,pNetwork->tau1,pNetwork->tau2,pNetwork->g);
	fopen_s(&pOutputFrequency,outfilename,"w");
	for(int i=0;i<neuron_num;++i)
		fprintf(pOutputFrequency,"%lf ",freq[i]);
	fprintf(pOutputFrequency,"\n");
	fprintf(pOutputFrequency,"The network frequency is %lf\n",networkfreq);

	delete [] freq;
}


double CalcSumCompExp(const double*_phase,const int _S){
	double re=0.0,im=0.0;
	for(int i=0;i<_S;++i){
		re+=cos(_phase[i]);
		im+=sin(_phase[i]);
	}
	return sqrt(pow(re,2.0)+pow(im,2.0))/_S;
}


void MeanPhase(const pInfoNetwork  pNetwork,const int _m,const int _n){
	//Create directory...
	char dir[50];
	sprintf_s(dir,"Connor%dHH%dML1_%dML2_%d",pNetwork->Connor_num,pNetwork->HH_num,pNetwork->ML1_num,pNetwork->ML2_num);
	_mkdir(dir);

	//Create file pointers...
	char Infilename[80];
	sprintf_s(Infilename,".\\%s\\tau1=%.3lftau2=%.3lfg=%.3lfSpikeTime.dat",dir,pNetwork->tau1,pNetwork->tau2,pNetwork->g);
	FILE *pInputSpikeTime;
	if(0!=fopen_s(&pInputSpikeTime,Infilename,"r")){
		OutputSpikeTime(pNetwork);
		fopen_s(&pInputSpikeTime,Infilename,"w");
	}

	
}

/*
void SyncBursting(const double&_I){
	//Create directory...
	char dir[50];
	sprintf_s(dir,"I=%.3lfmV",_I);
	_mkdir(dir);

	//Create file pointers...
	char Infilename[80];
	sprintf_s(Infilename,".\\I=%.3lf Tau1=%.3lf Tau2=%.3lf G=%.3lf\\SpikeTime_0.dat",_I,TAU1,TAU2,G);
	FILE *test;
	if(0!=fopen_s(&test,Infilename,"r")){
		OutputSpikeTime(_I);
	}

	FILE *InputSpikeTime[N];
	for(int i=0;i<N;++i){
		sprintf_s(Infilename,".\\I=%.3lf Tau1=%.3lf Tau2=%.3lf G=%.3lf\\SpikeTime_%d.dat",_I,TAU1,TAU2,G,i);
		fopen_s(&InputSpikeTime[i],Infilename,"w");
	}

}*/

/*Calculate amplitude of average vector field*/
double AverVecField(const pInfoNetwork pNetwork){
	/*Create directory...*/
	char dir[50];
	sprintf_s(dir,"Connor%dHH%dML1_%dML2_%d",pNetwork->Connor_num,pNetwork->HH_num,pNetwork->ML1_num,pNetwork->ML2_num);
	_mkdir(dir);

	/*Create file pointer...*/
	FILE *pInputTimeSeries;
	char Infilename[80];
	sprintf_s(Infilename,".\\%s\\tau1=%.3lftau2=%.3lfg=%.3lfTimeSeries.dat",dir,pNetwork->tau1,pNetwork->tau2,pNetwork->g);
	if(0!=fopen_s(&pInputTimeSeries,Infilename,"r")){
		printf("Can't open file in function AverVecField!");
		exit(1);
	}
	int neuron_num=pNetwork->Connor_num+pNetwork->HH_num+pNetwork->ML1_num+pNetwork->ML2_num;
	double t;
	double voltage,sum_voltage;
	double averVecFieldMin,averVecFieldMax;
	sum_voltage=0;
	if(fscanf_s(pInputTimeSeries,"%lf",&t)==EOF){
		printf("can't scanf !");
		exit(1);
	}
	for(int i=0;i<neuron_num;++i){
		fscanf_s(pInputTimeSeries,"%lf",&voltage);
		sum_voltage+=voltage;
	}
	averVecFieldMin=sum_voltage/neuron_num;
	averVecFieldMax=averVecFieldMin;
	while(fscanf_s(pInputTimeSeries,"%lf",&t)==EOF){
		sum_voltage=0;
		for(int i=0;i<neuron_num;++i){
			fscanf_s(pInputTimeSeries,"%lf",&voltage);
			sum_voltage+=voltage;
		}
		sum_voltage/=neuron_num;
		if(averVecFieldMax<sum_voltage) 
			averVecFieldMax=sum_voltage;
		if(averVecFieldMin>sum_voltage)
			averVecFieldMin=sum_voltage;
	}

	return (averVecFieldMax-averVecFieldMin);

}

/*Output tau1 tau2 average vector field*/
void Tau1Tau2AverVecField(int _Connor_num,int _HH_num,int _ML1_num, int _ML2_num){
	infoNetwork theinfoNetwork;
	theinfoNetwork.Connor_num=_Connor_num;
	theinfoNetwork.HH_num=_HH_num;
	theinfoNetwork.ML1_num=_ML1_num;
	theinfoNetwork.ML2_num=_ML2_num;
	pInfoNetwork ptheinfoNetwork=&theinfoNetwork;

	//Create directory...
	char dir[40];
	sprintf_s(dir,"Connor%dHH%dML1_%dML2_%d",_Connor_num,_HH_num,_ML1_num,_ML2_num);
	_mkdir(dir);
	//Create file pointer...
	FILE *pOutputAverVecField;
	char Outfilename[60];
	sprintf_s(Outfilename,".\\%s\\tau1_tau2_aver_vec_field.dat",dir);
	fopen_s(&pOutputAverVecField,Outfilename,"w");

	theinfoNetwork.g=0.1;
	double averVecField;
	for(double i=2;i<11;i+=1){
		for(double j=1;j<i;j+=1){
			theinfoNetwork.tau1=i;
			theinfoNetwork.tau2=j;
			averVecField=AverVecField(ptheinfoNetwork);
			fprintf(pOutputAverVecField,"%lf %lf %lf\n",i,j,averVecField);
		}
	}
}

/*Calculate synchronization factor*/
double SyncFactor(const pInfoNetwork pNetwork){
	/*Create directory...*/
	char dir[50];
	sprintf_s(dir,"Connor%dHH%dML1_%dML2_%d",pNetwork->Connor_num,pNetwork->HH_num,pNetwork->ML1_num,pNetwork->ML2_num);
	_mkdir(dir);

	/*Create file pointer...*/
	FILE *pInputTimeSeries;
	char Infilename[80];
	sprintf_s(Infilename,".\\%s\\tau1=%.3lftau2=%.3lfg=%.3lfTimeSeries.dat",dir,pNetwork->tau1,pNetwork->tau2,pNetwork->g);
	fopen_s(&pInputTimeSeries,Infilename,"r");
	double t,voltage;
	double sum_averVkSquare=0.0,sum_averVk=0.0;
	int nTime=0;
	int neuron_num=pNetwork->Connor_num+pNetwork->HH_num+pNetwork->ML1_num+pNetwork->ML2_num;
	double *sum_VikSquare=new double[neuron_num];
	double *sum_Vik=new double[neuron_num];
	memset(sum_VikSquare,0,sizeof(double)*neuron_num);
	memset(sum_Vik,0,sizeof(double)*neuron_num);
	double sumVecField,averVecField;
	while(fscanf_s(pInputTimeSeries,"%lf",&t)!=EOF){
		nTime++;
		sumVecField=0;
		for(int i=0;i<neuron_num;++i){
			fscanf_s(pInputTimeSeries,"%lf",&voltage);
			sumVecField+=voltage;
			sum_Vik[i]+=voltage;
			sum_VikSquare[i]+=(voltage*voltage);
		}
		averVecField=sumVecField/neuron_num;
		sum_averVk+=(averVecField);
		sum_averVkSquare+=(averVecField*averVecField);
	}
	fclose(pInputTimeSeries);
	double R,numerator=0,denominator=0;
	numerator=((sum_averVkSquare/nTime)-(sum_averVk/nTime)*(sum_averVk/nTime));
	for(int i=0;i<neuron_num;++i){
		denominator+=((sum_VikSquare[i]/nTime)-(sum_Vik[i]/nTime)*(sum_Vik[i]/nTime));
	}
	denominator/=neuron_num;

	R=numerator/denominator;

	delete []sum_Vik;
	delete []sum_VikSquare;

	return R;


}

/*Output tau1 tau2 synchronization factor*/
void Tau1Tau2SyncFactor(const int _Connor_num,const int _HH_num,const int _ML1_num,const int _ML2_num){
	infoNetwork theinfoNetwork;
	theinfoNetwork.Connor_num=_Connor_num;
	theinfoNetwork.HH_num=_HH_num;
	theinfoNetwork.ML1_num=_ML1_num;
	theinfoNetwork.ML2_num=_ML2_num;
	pInfoNetwork ptheinfoNetwork=&theinfoNetwork;

	//Create directory...
	char dir[40];
	sprintf_s(dir,"Connor%dHH%dML1_%dML2_%d",_Connor_num,_HH_num,_ML1_num,_ML2_num);
	_mkdir(dir);
	//Create file pointer...
	FILE *pOutputSyncFactor;
	char Outfilename[60];
	sprintf_s(Outfilename,".\\%s\\tau1_tau2_sync_factor.dat",dir);
	fopen_s(&pOutputSyncFactor,Outfilename,"w");

	theinfoNetwork.g=0.1;
	double R;
	for(double i=2;i<11;i+=0.2){
		for(double j=1;j<i;j+=0.2){
			theinfoNetwork.tau1=i;
			theinfoNetwork.tau2=j;
			R=SyncFactor(&theinfoNetwork);
			fprintf(pOutputSyncFactor,"%lf %lf %lf\n",i,j,R);
			printf("i=%lf j=%lf is completed!\n",i,j);
		}
	}
	fclose(pOutputSyncFactor);
}