#include "stdafx.h"

#include "ConnorModel.h"

const double V_NA=55.0;
const double V_K=-72.0;
const double V_L=-17.0;
const double V_A=-75.0;
const double G_NA=120.0;
const double G_K=20.0;
const double G_L=0.3;
const double G_A=47.7;
const double C=1.0;

Connor::Connor():No(0){
	strcpy_s(neuron_model,"ConnorModel");
}

Connor::Connor(const int&_No,bool _IsExcitation):No(_No){
	I=13; //just above Bifurcation point 
	strcpy_s(neuron_model,"ConnorModel");
}

Connor::~Connor(){}

void Connor::InitVar(){
	V=-40.0;
	m=0.0;
	h=0.0;
	n=0.0;
	A=0.0;
	B=0.0;
}

void Connor::RandVar(){
	V=-40+20.0*rand()/(double)RAND_MAX;
	m=0.0;
	h=0.0;
	n=0.0;
	A=0.0;
	B=0.0;
}

void Connor::equation_s(){
	DV=(I-G_NA*m*m*m*h*(V-V_NA)-G_K*n*n*n*n*(V-V_K)-G_L*(V-V_L)-G_A*A*A*A*B*(V-V_A))/C;
	Dm=alpha_m()-(alpha_m()+beta_m())*m;
	Dh=alpha_h()-(alpha_h()+beta_h())*h;
	Dn=alpha_n()-(alpha_n()+beta_n())*n;
	DA=(a_infinity()-A)/tau_a();
	DB=(b_infinity()-B)/tau_b();
}

void Connor::equation(){
	DV=(I-G_NA*m*m*m*h*(V-V_NA)-G_K*n*n*n*n*(V-V_K)-G_L*(V-V_L)-G_A*A*A*A*B*(V-V_A)+Network::sm_gCouple[No])/C;
	Dm=alpha_m()-(alpha_m()+beta_m())*m;
	Dh=alpha_h()-(alpha_h()+beta_h())*h;
	Dn=alpha_n()-(alpha_n()+beta_n())*n;
	DA=(a_infinity()-A)/tau_a();
	DB=(b_infinity()-B)/tau_b();
}

void Connor::euler_s(){
	equation_s();
	V+=DV*single_dt;
	m+=Dm*single_dt;
	h+=Dh*single_dt;
	n+=Dn*single_dt;
	A+=DA*single_dt;
	B+=DB*single_dt;
}

void Connor::euler(){
	equation();
	V+=DV*Network::sm_dt;
	m+=Dm*Network::sm_dt;
	h+=Dh*Network::sm_dt;
	n+=Dn*Network::sm_dt;
	A+=DA*Network::sm_dt;
	B+=DB*Network::sm_dt;
}

void Connor::rungekutta_s(){
	V0=V;m0=m;h0=h;n0=n;A0=A;B0=B;
	/*1.....*/
	equation_s();
	DV1=DV*single_dt;Dm1=Dm*single_dt;Dh1=Dh*single_dt;Dn1=Dn*single_dt;DA1=DA*single_dt;DB1=DB*single_dt;
	V=V0+DV1/2.0;
	m=m0+Dm1/2.0;
	h=h0+Dh1/2.0;
	n=n0+Dn1/2.0;
	A=A0+DA1/2.0;
	B=B0+DB1/2.0;
	/*2.....*/
	equation_s();
	DV2=DV*single_dt;Dm2=Dm*single_dt;Dh2=Dh*single_dt;Dn2=Dn*single_dt;DA2=DA*single_dt;DB2=DB*single_dt;
	V=V0+DV2/2.0;
	m=m0+Dm2/2.0;
	h=h0+Dh2/2.0;
	n=n0+Dn2/2.0;
	A=A0+DA2/2.0;
	B=B0+DB2/2.0;
	/*3.....*/
	equation_s();
	DV3=DV*single_dt;Dm3=Dm*single_dt;Dh3=Dh*single_dt;Dn3=Dn*single_dt;DA3=DA*single_dt;DB3=DB*single_dt;
	V=V0+DV3;
	m=m0+Dm3;
	h=h0+Dh3;
	n=n0+Dn3;
	A=A0+DA3;
	B=B0+DB3;
	/*4.....*/
	equation_s();
	DV4=DV*single_dt;Dm4=Dm*single_dt;Dh4=Dh*single_dt;Dn4=Dn*single_dt;DA4=DA*single_dt;DB4=DB*single_dt;
	V=V0+(DV1+2.0*DV2+2.0*DV3+DV4)/6.0;
	m=m0+(Dm1+2.0*Dm2+2.0*Dm3+Dm4)/6.0;
	h=h0+(Dh1+2.0*Dh2+2.0*Dh3+Dh4)/6.0;
	n=n0+(Dn1+2.0*Dn2+2.0*Dn3+Dn4)/6.0;
	A=A0+(DA1+2.0*DA2+2.0*DA3+DA4)/6.0;
	B=B0+(DB1+2.0*DB2+2.0*DB3+DB4)/6.0;
	/*Update time...*/
	//t+=sm_dt;
}

void Connor::rungekutta(){
	V0=V;m0=m;h0=h;n0=n;A0=A;B0=B;
	/*1.....*/
	equation();
	DV1=DV*Network::sm_dt;Dm1=Dm*Network::sm_dt;Dh1=Dh*Network::sm_dt;Dn1=Dn*Network::sm_dt;DA1=DA*Network::sm_dt;DB1=DB*Network::sm_dt;
	V=V0+DV1/2.0;
	m=m0+Dm1/2.0;
	h=h0+Dh1/2.0;
	n=n0+Dn1/2.0;
	A=A0+DA1/2.0;
	B=B0+DB1/2.0;
	/*2.....*/
	equation();
	DV2=DV*Network::sm_dt;Dm2=Dm*Network::sm_dt;Dh2=Dh*Network::sm_dt;Dn2=Dn*Network::sm_dt;DA2=DA*Network::sm_dt;DB2=DB*Network::sm_dt;
	V=V0+DV2/2.0;
	m=m0+Dm2/2.0;
	h=h0+Dh2/2.0;
	n=n0+Dn2/2.0;
	A=A0+DA2/2.0;
	B=B0+DB2/2.0;
	/*3.....*/
	equation();
	DV3=DV*Network::sm_dt;Dm3=Dm*Network::sm_dt;Dh3=Dh*Network::sm_dt;Dn3=Dn*Network::sm_dt;DA3=DA*Network::sm_dt;DB3=DB*Network::sm_dt;
	V=V0+DV3;
	m=m0+Dm3;
	h=h0+Dh3;
	n=n0+Dn3;
	A=A0+DA3;
	B=B0+DB3;
	/*4.....*/
	equation();
	DV4=DV*Network::sm_dt;Dm4=Dm*Network::sm_dt;Dh4=Dh*Network::sm_dt;Dn4=Dn*Network::sm_dt;DA4=DA*Network::sm_dt;DB4=DB*Network::sm_dt;
	V=V0+(DV1+2.0*DV2+2.0*DV3+DV4)/6.0;
	m=m0+(Dm1+2.0*Dm2+2.0*Dm3+Dm4)/6.0;
	h=h0+(Dh1+2.0*Dh2+2.0*Dh3+Dh4)/6.0;
	n=n0+(Dn1+2.0*Dn2+2.0*Dn3+Dn4)/6.0;
	A=A0+(DA1+2.0*DA2+2.0*DA3+DA4)/6.0;
	B=B0+(DB1+2.0*DB2+2.0*DB3+DB4)/6.0;
	/*Update time...*/
	//t+=sm_dt;
}