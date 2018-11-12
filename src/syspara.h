//#ifndef __SYSPARA_H_INCLUDE 
//#define __SYSPARA_H_INCLUDE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "mkl.h"
#include "/home/tsumoto/lib/xhplot.h"

#define NN 22
#define BUF 200
#define NUM 20

#define R 8314.472
#define F 96485.33771638995
#define T 310

#define dvm 5
#define Emax 2000
#define Emin -2000
#define VNMAX (Emax-Emin)*dvm+1

struct varstruct {

	int datas;
	int line_wid[NUM];
	
	int n;
	double Istim,dIstim;
	double coef,dist;
	double s1,s2;
	double RTonF,RTon2F;

	// Drug concentration for Nifekalant
	double drug_conc,ic50,ec50,hill,ehill;
	double normal_block,fraction_facil;
	double block_rate;

	// Model & Cell type & Drug type
	// Cell tupe
	// 0: Endo, 1: Mid, 2: Epi
	int model_type,celltype,drug_type;

	// Cell Geometry
	double RGC;
	double length,a;
	double vcell,ageo;
	double acap,vmyo;
	double vmito,vsr;
	double vnsr,vjsr;
	double vcleft;
	double vr1,vr2,vr3,vr4,vr5;

	// Ion Valences 
	double zna,zk,zca;

	// time ratio
	double ndis;

	// Ion concentrations
	//	double nao,ko,cao;

	// Tablization of exp functions
	// Fast sodium current
	double *Tam,*Tbm,*Tah,*Tbh,*Taj,*Tbj,*TEna;
	
	// L-type calcium current
	double *TexpCa,*TexpNa,*TexpK;
	double *Tdss,*Ttaud,*Tfss,*Ttauf;

	// T-type calcium current
	double *Tbss,*Ttaub,*Tgss,*Ttaug;

	// Rapidly activating potassium current
	double *Txrss,*Ttauxr,*Tr,*TEk,*Txr2ss;
	
	// Slowly activating potassium current
	double *Txs1ss,*Ttauxs1;
	
	// Inward rectificating potassium current
	//double *T1Ki,*T2Ki,*T3Ki,*T4Ki;
	//double *T5Ki,*T6Ki,*T7Ki,*T8Ki;

	// Plateau Potassium current
	double *Tkp;
	
	// / Na-Activated K Channel 
	double *Tpov;

	// Transient outward current
	double *Trvdv,*Tazdv,*Tbzdv,*Taydv,*Tbydv;

	// using Inak and Inaca
	double *Texp0,*Texp1,*Texp2,*Texp3,*Tnak;

	// using Irel and Irel_ol
	double *Tc,*Tol;

	// Change rates for K channel conductances
	double ikrf,iksf,ikxf,itof;
	
	// reversal potential
	double Ena,Ek,Eca,Eks;

	// Fast sodium current
	double ina,gna,am,bm,ah,bh,aj,bj;

	// L-type calcium current
	double ilcatot;
	double ibarca,ibarna,ibark;
	double ilca,ilcana,ilcak;
	double dss,taud;                  // Steady-state value of activation gate d and time const
	double fss,tauf;     // Steady-state value of inactivation gate f and time const
	double kmca;
	double pca,gacai,gacao;
	double pna,ganai,ganao;           // Permiability of membrane to Na (cm/s)
	double pk,gaki,gako;              // Permiability of membrane to K (cm/s)
	double fca;
	double gca_rate;

	// T-type calcium current
	double icat,gcat,bss,taub,gss,taug;

	// Rapidly activating potassium current
	double gkr,xrss,tauxr,r;
	double xr2ss;
	double ikr,ikr_with,ikr_without;
	double gkr_max;
	
	// Slowly activating potassium current
	double iks,gks,xs1ss,tauxs1,xs2ss,tauxs2;
	double prnak;
	double gks_rate;
	double gks_max;

	// Potassium current IK1 (time-dependent)
	double iki,gki,aki,bki,kin;
	double gk1_rate;

	// Plateau Potassium current
	double ikp,gkp,kp;
	double gkp_max;

	// Na-activated potassium channel
		double ikna,gkna,pona,pov,nkna,kdkna;

	// ATP-sensitive potassium channel
	double gkatp,gkbaratp,patp,natp;
	double nicholsarea,atpi,hatp,katp;
	double ikatp;

	// Ito Transient Outward Current
	// (Dumaine et al. Circ Res 1999;85:803-809)
	double gitodv,gitodv_max,rvdv;
	double azdv,bzdv,aydv,bydv;
	double Ekdv,gto,ito;

	// Sodium-Calcium Exchanger V-S
	double c1,c2,gammas,inaca;

	// Sodium-Potassium Pump
	double fnak,sigma,ibarnak,kmnai,kmko;
	double inak;
	
	// Nonspecific Ca-activated Current 
	double ibarnsna,ibarnsk,pnsca,kmnsca;
	double insna,insk;

	// Sarcolemmal Ca Pump
	double ibarpca,kmpca,ipca;

	// Ca Background Current
	double icab,gcab;

	// Na Background Current
	double inab,gnab;

	// Ultra rapid Potassium current
	double gkur,ikur;

	// Total Ion currents 
	double INa_total;
	double IK_total;
	double ICa_total;
	
	// Difference total Ion current 
	double dICa_total;
	
	// NSR Ca Ion Concentration Changes
	double iup,ileak;      // Ca uptake from myo. to NSR (mM/ms)
	double kmup,iupbar,nsrbar;

	// JSR Ca Ion Concentration Changes
	int boolien;
	double Irel_cicr,Irel_jsr_overload;
	double dICa_total_new,dCa_ion,ICa_total_old;
	double gmaxrel;
	double jsr_new,grelbarjsrol,tauon,tauoff;
	double csqnbar,kmcsqn,t_cicr,t_overload;
	double bjsr,cjsr,csqn,csqnth,swspontan;
	double djsr,jsr;

	// Translocation of Ca Ions from NSR to JSR
	double tautr,itr;

	// Ca concentration
	double cmdnbar,trpnbar,kmtrpn,kmcmdn;
	double trpn,cmdn;
	double Ca_total,gpig;
	double bmyo,cmyo,dmyo;

	// Extracellular ion concentrations
	double Na_bulk,K_bulk,Ca_bulk;
	double tau_diff;

	// Base Currnt Stimulus
	double Istim_base;

	// Sttimulus parameters
	double BCL;  // Base cycle length = stimulus period
	double dt;  // time step
	int beat; // Number of stimulus

	// debug variable
	double ca_pre,dca_now;

    int m;
    int l;

	double x0[NUM][NN];
    double tsign[NUM];
    double tend[NUM];

	int pflag;
	int write, graph;
	int write0;
	int half;
	int deb;
	int pswitch, sswitch;
	int out_data;

} var;

// for main
void make_ExPTable();
void eular(int n, double h, double x[],double t);
void function(double x[],double f[],double t);
//void input_para(FILE *);
void static_paras(FILE *);
void mem();
void close_mem();

void eventloop(FILE *, int *mode, int *P, double m[]);
void orbit(int *mode, double m[], double x2);
void draw_p(int *mode, int P, double x[], double x2);
void mouse(int *mode, double x[], double x2);

// for dataout
void data_out(FILE *, double t, double u[]);
void current(FILE *,FILE *,FILE *, FILE *, double t,double x[]);

void out_ikr (FILE *, double time, double p[]);
void out_iks (FILE *, double time, double p[]);
void out_ical(FILE *, double time, double p[]);
void out_inaca (FILE *, double time, double p[]);
void out_inak (FILE *, double time, double p[]);
void out_cicr (FILE *, double time, double p[]);

// for calculation of Ion currents
void comp_rev(double x[]);
void comp_ina(double x[]);
void comp_icat(double x[]);
void comp_ical(double x[]);
void comp_ikr(double x[]);
void comp_iki(double x[]);
void comp_iks(double x[]);
void comp_ikp(double x[]);
void comp_ikur(double x[]);
void comp_ikna(double x[]);
void comp_ikatp(double x[]);
void comp_ito(double x[]);
void comp_inaca(double x[]);
void comp_inak(double x[]);
void comp_insca(double x[]);
void comp_ipca(double x[]);
void comp_icab(double x[]);
void comp_inab(double x[]);
void comp_iconcent (double x[]);
void comp_iconcent2 (double x[]);
void conc_nsr(double x[]);
void conc_jsr(double x[]);
void conc_itr(double x[]);
void conc_cai(double x[]);
void conc_cleft(double x[]);

//void current_ikr(FILE *, double t, double x[]);
//void current_iks(FILE *, double t, double x[]);
//void current_ical(FILE *, double t, double x[]);
//void current_incx(FILE *, double t, double x[]);
//void current_inak(FILE *, double t, double x[]);
//void current_ik1(FILE *, double t, double x[]);
//void current_ina(FILE *, double t, double x[]);
//void current_ito(FILE *, double t, double x[]);
//void current_it(FILE *, double t, double x[]);
//void current_irel(FILE *, double t, double x[]);

//void main(int argc, char **argv);

