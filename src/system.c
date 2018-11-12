
#include "syspara.h"

inline void function(double x[], double f[], double t)
{
	
	comp_rev(x);
	comp_ina(x);
	comp_ical(x);
	comp_icat(x);
	comp_ikr(x);
	comp_iks(x);
	comp_iki(x);
	comp_ikp(x);
	comp_ikna(x);
	comp_ikatp(x);
	comp_ito(x);
	comp_inaca(x);
	comp_inak(x);
	comp_insca(x);
	comp_ipca(x);
	comp_icab(x);
	comp_inab(x);

	var.INa_total = var.ina + var.inab + var.ilcana + var.insna + 3.0*var.inak + 3.0*var.inaca;
	var.IK_total = var.ikr + var.iks + var.iki + var.ikp + var.ilcak + var.insk - 2.0*var.inak + var.ito + var.ikna + var.ikatp + var.Istim;
	var.ICa_total = var.ilca + var.icab + var.ipca - 2.0*var.inaca + var.icat;
	
	comp_iconcent(x);

	f[0] = -(var.INa_total+var.IK_total+var.ICa_total);
	f[1] = var.am*(1.0 - x[1]) - var.bm*x[1]; // m
	f[2] = var.ah*(1.0 - x[2]) - var.bh*x[2]; // h
	f[3] = var.aj*(1.0 - x[3]) - var.bj*x[3]; // j
	f[4] = (var.dss - x[4])/var.taud;	// dss
	f[5] = (var.fss - x[5])/var.tauf;	// fss

	f[6] = (var.bss - x[6])/var.taub;	//bss
	f[7] = (var.gss - x[7])/var.taug;	//gss

	f[8] = (var.xrss - x[8])/var.tauxr;	//xrss

	f[9] = (var.xs1ss - x[9])/var.tauxs1;	//xs1ss
	f[10] = (var.xs2ss - x[10])/var.tauxs2;	//xs2ss

	f[11] = var.azdv*(1.0 - x[11]) - var.bzdv*x[11];
	f[12] = var.aydv*(1.0 - x[12]) - var.bydv*x[12];

	f[13] = -var.INa_total*(var.acap/(var.vmyo*var.zna*F));
	f[14] = -var.IK_total*(var.acap/(var.vmyo*var.zk*F));
	f[15] = -var.ICa_total*(var.acap/(var.vmyo*var.zca*F)) - ((var.iup-var.ileak)*var.vnsr/var.vmyo) + (var.Irel_cicr+var.Irel_jsr_overload)*(var.vjsr/var.vmyo);

	f[16] = var.iup - var.ileak - var.itr*var.vjsr/var.vnsr;
	f[17] = var.itr - var.Irel_cicr - var.Irel_jsr_overload;

	f[18] = (var.Na_bulk - x[18])/var.tau_diff + var.INa_total*var.acap/(var.vcleft*F);
	f[19] = (var.K_bulk - x[19])/var.tau_diff + var.IK_total*var.acap/(var.vcleft*F);
	f[20] = (var.Ca_bulk - x[20])/var.tau_diff + var.ICa_total*var.acap/(var.vcleft*var.zca*F);

	f[21] = (var.xr2ss - x[21])/var.tauxr;       //xrss

}

void comp_rev(double x[])
{
	var.Ena = var.RTonF*log(x[18]/x[13]);
	var.Ek = var.RTonF*log(x[19]/x[14]);
	var.Eks = ((R*T)/F)*log((x[19]+var.prnak*x[18])/(x[14]+var.prnak*x[15]));
	var.Eca = (var.RTonF/var.zca)*log(x[20]/x[15]); // [Ca]i = x[15],[Ca]o = x[20]

}

// Fast Sodium Current (time dependant) */
inline void comp_ina(double x[])
{

	MKL_INT iV=0, iNa=0;
	double V1,V2,d1,d2;
	double Na1,Na2,k1,k2;
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.am = var.Tam[iV]*d2 + var.Tam[iV+1]*d1;
	var.bm = var.Tbm[iV]*d2 + var.Tbm[iV+1]*d1;
	var.ah = var.Tah[iV]*d2 + var.Tah[iV+1]*d1;
	var.bh = var.Tbh[iV]*d2 + var.Tbh[iV+1]*d1;
	var.aj = var.Taj[iV]*d2 + var.Taj[iV+1]*d1;
	var.bj = var.Tbj[iV]*d2 + var.Tbj[iV+1]*d1;

	var.ina = var.gna*x[1]*x[1]*x[1]*x[2]*x[3]*(x[0]-var.Ena);

}

// L-type calcium current
inline void comp_ical(double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	double expCa,expNa,expK;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.dss = var.Tdss[iV]*d2 + var.Tdss[iV+1]*d1;
	var.taud = var.Ttaud[iV]*d2 + var.Ttaud[iV+1]*d1;
	var.fss= var.Tfss[iV]*d2 + var.Tfss[iV+1]*d1;
	var.tauf= var.Ttauf[iV]*d2 + var.Ttauf[iV+1]*d1;

	expCa = var.TexpCa[iV]*d2 + var.TexpCa[iV+1]*d1;
	expNa = var.TexpNa[iV]*d2 + var.TexpNa[iV+1]*d1;
	expK = var.TexpK[iV]*d2 + var.TexpK[iV+1]*d1;

	var.ibarca = var.pca*var.zca*var.zca*((x[0]*F*F)/(R*T))*((var.gacai*x[15]*expCa-var.gacao*x[20])/(expCa-1.0));
	var.ibarna = var.pna*var.zna*var.zna*((x[0]*F*F)/(R*T))*((var.ganai*x[13]*expNa-var.ganao*x[18])/(expNa-1.0));
	var.ibark = var.pk*var.zk*var.zk*((x[0]*F*F)/(R*T))*((var.gaki*x[14]*expK-var.gako*x[19])/(expK-1.0));

	var.fca=1.0/(1.0+x[15]/var.kmca);

	var.ilca = x[4]*x[5]*var.fca*var.ibarca;
	var.ilcana = x[4]*x[5]*var.fca*var.ibarna;
	var.ilcak = x[4]*x[5]*var.fca*var.ibark;

	var.ilcatot = var.ilca + var.ilcana + var.ilcak;

}

// Current through T-type Ca Channel */
inline void comp_icat(double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.bss = var.Tbss[iV]*d2 + var.Tbss[iV+1]*d1;
	var.taub = var.Ttaub[iV]*d2 + var.Ttaub[iV+1]*d1;
	var.gss = var.Tgss[iV]*d2 + var.Tgss[iV+1]*d1;
	var.taug = var.Ttaug[iV]*d2 + var.Ttaug[iV+1]*d1;
	
	var.icat = var.gcat*x[6]*x[6]*x[7]*(x[0]-var.Eca);

}

// Rapidly Activating Potassium Current 
inline void comp_ikr (double x[])
{

	
	MKL_INT iV=0,iK=0;   
	double V1,V2,d1,d2;
	double K1,K2,dd1,dd2;

	double Ikr,Ikr2;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	// normal (unfacilitated fraction)
	var.xrss = var.Txrss[iV]*d2 + var.Txrss[iV+1]*d1;
	var.tauxr = var.Ttauxr[iV]*d2 + var.Ttauxr[iV+1]*d1;
	var.r = var.Tr[iV]*d2 + var.Tr[iV+1]*d1;

	// facilitated fraction
	var.xr2ss = var.Txr2ss[iV]*d2 + var.Txr2ss[iV+1]*d1;

	var.gkr = 0.02614*sqrt(x[19]/5.4);
	
	Ikr  = var.gkr*x[8]*var.r*(x[0]-var.Ek);
	Ikr2 = var.gkr*x[21]*var.r*(x[0]-var.Ek);

	// Drug effects
	if(var.drug_type == 0){ // without facilitation effect

		var.ikr_without = var.normal_block*Ikr;
		var.ikr_with    = 0.0;

	} else if(var.drug_type == 1){ // with facilitation effect

		var.ikr_without = var.normal_block*var.fraction_facil*Ikr;
		var.ikr_with    = var.normal_block*(1.0-var.fraction_facil)*Ikr2;
	}
	
	var.ikr = var.ikr_with + var.ikr_without;

}

// Slowly Activating Potassium Current 
inline void comp_iks (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	double gks;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.xs1ss = var.Txs1ss[iV]*d2 + var.Txs1ss[iV+1]*d1;
	var.xs2ss = var.xs1ss;
	var.tauxs1 = var.Ttauxs1[iV]*d2 + var.Ttauxs1[iV+1]*d1;
	var.tauxs2 = 4.0*var.tauxs1;

	gks = var.iksf*(1.0+0.6/(1.0+pow((0.000038/x[15]),1.4)));

	var.iks = var.gks_rate*gks*x[9]*x[10]*(x[0]-var.Eks);

}

// Potassium Current (time-independant) Ik1
inline void comp_iki (double x[])
{
	
	var.aki = 1.02/(1.0+exp(0.2385*(x[0]-var.Ek-59.215)));
	var.bki = (0.49124*exp(0.08032*(x[0]-var.Ek+5.476))+exp(0.06175*(x[0]-var.Ek-594.31)))/(1.0+exp(-0.5143*(x[0]-var.Ek+4.753)));

	var.kin = var.aki/(var.aki+var.bki);

	//var.gki = var.gk1_rate*0.75*sqrt(x[22]/5.4);
	var.gki = 0.75*sqrt(x[19]/5.4);

	var.iki = var.gki*var.kin*(x[0]-var.Ek);

}

// Plateau Potassium Current
inline void comp_ikp (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.kp = var.Tkp[iV]*d2 + var.Tkp[iV+1]*d1;
	
	var.ikp = var.gkp*var.kp*(x[0]-var.Ek);

}

// Na-Activated K Channel 
inline void comp_ikna (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.pov = var.Tpov[iV]*d2 + var.Tpov[iV+1]*d1;
	var.pona = 0.85/(1.0 + pow((var.kdkna/x[13]),2.8));
	
	//var.ikna = var.gkna*var.pona*var.pov*(x[0]-var.Ek);
	var.ikna = 0.0;
}

// ATP-Sensitive K Channel
// Note: If you wish to use this current in your simulations,
// there are additional changes which must be made to the code 
// as detailed in Cardiovasc Res 1997;35:256-272 
// [K]o = x[19] extracellular K ion concentration
inline void comp_ikatp (double x[])
{

	var.gkbaratp = var.gkatp*var.patp*(pow((x[19]/4.0),var.natp));
	//var.ikatp = var.gkbaratp*(x[0]-var.Ek);
	var.ikatp = 0.0;

}


// Ito Transient Outward Current
// Dumaine et al. Circ Res 1999;85:803-809
inline void comp_ito (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.rvdv = var.Trvdv[iV]*d2 + var.Trvdv[iV+1]*d1;
	var.azdv = var.Tazdv[iV]*d2 + var.Tazdv[iV+1]*d1;
	var.bzdv = var.Tbzdv[iV]*d2 + var.Tbzdv[iV+1]*d1;
	var.aydv = var.Taydv[iV]*d2 + var.Taydv[iV+1]*d1;
	var.bydv = var.Tbydv[iV]*d2 + var.Tbydv[iV+1]*d1;

	//var.ito = var.gto*x[11]*x[11]*x[11]*x[12]*var.rvdv*(x[0]-var.Ek);
	var.ito = 0.0;

}

// Sodium-Calcium Exchanger V-S
inline void comp_inaca (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	double exp2,exp3; 

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	exp2=var.Texp2[iV]*d2 + var.Texp2[iV+1]*d1;
	exp3=var.Texp3[iV]*d2 + var.Texp3[iV+1]*d1;

    var.inaca = var.c1*exp2*((exp3*x[13]*x[13]*x[13]*x[20]-x[18]*x[18]*x[18]*x[15])/(1.0+var.c2*exp2*(exp3*x[13]*x[13]*x[13]*x[20]+x[18]*x[18]*x[18]*x[15])));

}

// Sodium-Potassium Pump NaK ATPase
inline void comp_inak (double x[])
{

	MKL_INT iV=0,iNa=0;
	double V1,V2,d1,d2;
	double Na1,Na2,k1,k2;
	double exp0,exp1,tab_nak; 

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	exp0=var.Texp0[iV]*d2 + var.Texp0[iV+1]*d1;
	exp1=var.Texp1[iV]*d2 + var.Texp1[iV+1]*d1;

	var.sigma = (exp(x[21]/67.3)-1.0)/7.0;
	var.fnak = 1.0/(exp0+0.0365*var.sigma*exp1);
	tab_nak = 1.0/(1.0+pow(var.kmnai/x[13],2.0));

	var.inak = var.ibarnak*var.fnak*tab_nak*(x[19]/(x[19]+var.kmko));

}

// Nonspecific Ca-activated Current
inline void comp_insca (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	double expNa,expK;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	expNa = var.TexpNa[iV]*d2 + var.TexpNa[iV+1]*d1;
	expK = var.TexpK[iV]*d2 + var.TexpK[iV+1]*d1;

	var.ibarnsna = var.pnsca*var.zna*var.zna*((x[0]*F*F)/(R*T))*((var.ganai*x[13]*expNa-var.ganao*x[18])/(expNa-1.0));
	var.ibarnsk = var.pnsca*var.zk*var.zk*((x[0]*F*F)/(R*T))*((var.gaki*x[14]*expK-var.gako*x[19])/(expK-1.0));

	//var.insna = var.ibarnsna/(1.0 + pow(var.kmnsca/x[15],3.0));
	//var.insk = var.ibarnsk/(1.0 + pow(var.kmnsca/x[15],3.0));
	var.insna = 0.0;
	var.insk = 0.0;

}

// Sarcolemmal Ca Pump 
inline void comp_ipca (double x[])
{
	var.ipca = (var.ibarpca*x[15])/(var.kmpca+x[15]);
}

// Ca Background Current 
inline void comp_icab (double x[])
{
	var.icab = var.gcab*(x[0]-var.Eca);
}

// Na Background Current 
inline void comp_inab (double x[])
{
	var.inab = var.gnab*(x[0] - var.Ena);
}

// X-L Du et al, J.Mol.Cell.Cardiol, 35,293-300, 2003
inline void comp_ikur (double x[])
{
    var.ikur = 0;
}

// The units of dnai is in mM.  Note that naiont should be multiplied by the
// cell capacitance to get the correct units.  Since cell capacitance = 1 uF/cm^2,
// it doesn't explicitly appear in the equation below.
// This holds true for the calculation of dki and dcai.

// Na concentration dose not use a function.
// K concentration dose not use a function.
inline void comp_iconcent (double x[])
{

	conc_nsr(x);
	conc_jsr(x);
	conc_itr(x);
	conc_cai(x);

}

inline void conc_nsr (double x[])
{

	double kleak;

	kleak = var.iupbar/var.nsrbar;
	var.ileak = kleak*x[16];

	var.iup = var.iupbar*x[15]/(x[15]+var.kmup);

}

inline void conc_jsr (double x[])
{

	double magrel;
	double on,off;
	double greljsrol;
	double fexp,f_ol;

	var.dICa_total_new = (var.ICa_total - var.ICa_total_old)/var.dt;

	if( x[0] > -35.0 ){
		if(var.dICa_total_new < var.dICa_total){
			//printf("reset;CICR2\n");
			if(var.boolien == 0){
				var.boolien = 1;
				var.t_cicr = 0.0;
			//	printf("reset;CICR\n");
			}
		}
	}

	fexp=exp((-var.t_cicr+4.0)/0.5);
	on = 1.0/(1.0+fexp);
	off = (1.0-1.0/(1.0+fexp));

	magrel = 1.0/(1.0+exp(((var.ICa_total)+5.0)/0.9));

	var.Irel_cicr = var.gmaxrel*on*off*magrel*(x[17]-x[15]);
	
	var.t_cicr += var.dt;

	//if(var.grelbarjsrol=!0){
	if(var.csqn>=var.csqnth){
		var.grelbarjsrol=4.0;
		var.t_overload = 0.0;
	} else {
		var.grelbarjsrol = 0.0;
		var.t_overload += var.dt;
	}
	
	//if(var.csqn>=var.csqnth){
	//	f_ol=(1.0-exp(-var.t_overload/var.tauon))*exp(-var.t_overload/var.tauoff);

	//	var.Irel_jsr_overload = var.grelbarjsrol*f_ol*(x[17]-x[15]);
	//} else {
		var.Irel_jsr_overload = 0;
//		var.t_overload += 0;
//	}

	var.csqn = var.csqnbar*(x[17]/(x[17] + var.kmcsqn));

}

inline void conc_itr (double x[])
{
	var.itr = (x[16]-x[17])/var.tautr; 
}

inline void conc_cai (double x[])
{

	var.trpn = var.trpnbar*(x[15]/(x[15]+var.kmtrpn));
	var.cmdn = var.cmdnbar*(x[15]/(x[15]+var.kmcmdn));

}

inline void comp_iconcent2 (double x[])
{
// Update Ca_JSR concentration x[20] = [Ca]_JSR
	var.bjsr = var.csqnbar - var.csqn - x[17] + var.kmcsqn;
	var.cjsr = var.kmcsqn*(var.csqn + x[17]);
	x[17] = (sqrt(var.bjsr*var.bjsr+4.0*var.cjsr)-var.bjsr)/2.0;

// Update intracellular Ca concentration x[18] = [Ca]_i

	var.Ca_total = var.trpn+var.cmdn+x[15];
	var.bmyo = var.cmdnbar+var.trpnbar+var.kmtrpn+var.kmcmdn - var.Ca_total;
	var.cmyo = (var.kmcmdn*var.kmtrpn)-(var.Ca_total*(var.kmtrpn+var.kmcmdn))+(var.trpnbar*var.kmcmdn)+(var.cmdnbar*var.kmtrpn);
	var.dmyo = -var.kmtrpn*var.kmcmdn*var.Ca_total;
	var.gpig = sqrt(var.bmyo*var.bmyo-3.0*var.cmyo);

	x[15] = (2.0*var.gpig/3.0)*cos(acos((9.0*var.bmyo*var.cmyo
			-2.0*var.bmyo*var.bmyo*var.bmyo-27.0*var.dmyo)/
			(2.0*pow((var.bmyo*var.bmyo-3.0*var.cmyo),1.5)))/3.0)-(var.bmyo/3.0);

}



