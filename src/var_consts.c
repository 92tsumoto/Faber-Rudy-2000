#include "syspara.h"

void static_paras(FILE *fp1)
{
	int i,j,k;

	// general
		var.RTonF = R*T/F;
		var.RTon2F = R*T/(2.0*F);

	// Cell Geometry */

		var.RGC = 2;  // The ratio of the actual surface area to the geometrical surface area (cm)
		var.length = 0.01;  // Length of the cell (cm)
		var.a = 0.0011;  // Length of the cell (cm)
		var.vcell = 1000*M_PI*var.a*var.a*var.length; // Cell Volume:3.801e-5 (uL)
		var.ageo = 2.0*M_PI*var.a*var.a+2.0*M_PI*var.a*var.length; // geometric membrane area: 7.671e-5 (cm^2)
		var.acap = var.ageo*var.RGC;	// Capacitive membrane area: 1.534e-4 cm^2 (cm^2)
		var.vmyo = var.vcell*0.68;      // Myoplasm volume (uL) = 68% for Cell volume
		var.vmito = var.vcell*0.26;     // Mitochondria volume (uL) = 26% for cell volume
		var.vsr = var.vcell*0.06;       // SR volume (uL)
		var.vnsr = var.vcell*0.0552;    // NSR volume (uL)
		var.vjsr = var.vcell*0.0048;    // JSR volume (uL)
		var.vcleft = var.vcell*0.12/0.88;  // Cleft volume (uL)

	// Ion Valences
		var.zna = 1;  // Na valence
		var.zk = 1;   // K valence
		var.zca = 2;  // Ca valence
	
	// Bulk ion concentrations
		var.Na_bulk = 140.0;     // Na (mM) Correct value
		//var.K_bulk = 4.5;      // K (mM)
		var.K_bulk = 5.4;      // K (mM)
		var.Ca_bulk = 1.8;     // Ca (mM)

	// classification by cell type
	//	var.gkr = 0.02614;                      // real value: gkr*ikrf
	//	var.gkr_max = var.gkr*var.ikrf;
	//	var.gks = 0.433;                        // real value: gks_max*iksf
	//	var.gks_max = var.gks*var.iksf; 
	//	var.gkp = 0.02;                         // real value: gkp_max*ikxf
	//	var.gkp_max = var.gkp*var.ikxf; 
	//	var.gitodv = 0.5;                       // real value: gitodv*itof
	//	var.gitodv_max = var.gitodv*var.itof; 
		
	if(var.celltype==0){ // Endo cell
		//var.ikrf = 0.6400; var.iksf = 0.1903; var.ikxf = 1.2699; var.itof = 0.1131;
		var.iksf = 0.289;
	} else if(var.celltype==1){	// Mid cell
		//var.ikrf = 0.5382; var.iksf = 0.0951; var.ikxf = 1.2699; var.itof = 0.1600;
		var.iksf = 0.125;
	} else if(var.celltype==2){ // Epi cell
		//var.ikrf = 1.2800; var.iksf = 0.1903; var.ikxf = 0.5040; var.itof = 0.2263;
		var.iksf = 0.433;
	}	

	// Sodium channel current
		var.gna = 16.0;	//(mS/uF);

	// L-type calcium current
		var.kmca = 0.0006;	// Half-saturation concentration of Ca channel (mM)
		//var.pca = 0.00054;	// Permiability of membrane to Ca (cm/s)
		var.pca = 0.5*0.00054;	// Permiability of membrane to Ca (cm/s)
		var.gacai = 1.0;	// Activity coefficient of Ca
		var.gacao = 0.341;	// Activity coefficient of Ca
		var.pna = 0.000000675;	// Permiability of membrane to Na (cm/s)
		var.ganai = 0.75;	// Activity coefficient of Na 
		var.ganao = 0.75;	// Activity coefficient of Na 
		var.pk = 0.000000193;  // Permiability of membrane to K (cm/s)
		var.gaki = 0.75;	// Activity coefficient of K 
		var.gako = 0.75;	// Activity coefficient of K 
	
	// T-type calcium channel current
		var.gcat = 0.05;

	// Rapid Activated Potassium Current: IKr
		//var.gkr = 0.02614*sqrt(var.ko/5.4);

	// Slow Activated Potassium Current: IKs
		var.prnak = 0.01833;     // Na/K Permiability Ratio

	// Inward rectifier potassium current: IK1
		//var.gki = 0.75*sqrt(var.ko/5.4);

	// Plateau Potassium Current
		var.gkp = 0.00552;

	// Na-Activated K Channel 
		var.gkna = 0.12848; // Maximum conductance (mS/uF) 
		var.nkna = 2.8; // Hill coefficient for Na dependance 
		var.kdkna = 66.0; // Dissociation constant for Na dependance(mM)

	// ATP-Sensitive K Channel
		var.natp = 0.24; // K dependence of ATP-sensitive K current 
		var.nicholsarea = 0.00005; // Nichol's area (cm^2) 
		var.atpi = 3.0; // Intracellular ATP concentraion (mM) 
		var.hatp = 2.0; // Hill coefficient 
		var.katp = 0.250; // Half-maximal saturation point of ATP-sensitive K current (mM)
		var.gkatp = 0.000195/var.nicholsarea;
		var.patp = 1.0/(1.0+(pow((var.atpi/var.katp),var.hatp)));

	// Ito Transient Outward Current
		var.gto = 0.5;		// (mS/uF); <-- LV case, else 1.1 (RV case)
		
	// Sodium-Calcium Exchanger V-S (NCX)
		var.c1 = 0.00025;   // Scaling factor for inaca (uA/uF)
		var.c2 = 0.0001;    // Half-saturation concentration of NaCa exhanger (mM)
		var.gammas = 0.15;  // Position of energy barrier controlling voltage dependance of inaca

	// Sodium-Potassium Pump (NaK ATPase)
		var.ibarnak = 2.25;   // Max. current through Na-K pump (uA/uF)
		var.kmnai = 10;    // Half-saturation concentration of NaK pump (mM)
		var.kmko = 1.5;    // Half-saturation concentration of NaK pump (mM)

	// Nonspecific Ca-activated Current
		var.pnsca = 0.000000175;   // Permiability of channel to Na and K (cm/s)
		var.kmnsca = 0.0012;       // Half-saturation concentration of NSCa channel (mM)
	
	// Sarcolemmal Ca Pump
		var.ibarpca = 1.15;  // Max. Ca current through sarcolemmal Ca pump (uA/uF)
		var.kmpca = 0.0005;  // Half-saturation concentration of sarcolemmal Ca pump (mM)

	// Ca Background Current 
		var.gcab = 0.003016; // Max. conductance of Ca background (mS/uF)

	// Na Background Current 
		var.gnab = 0.004;    // Max. conductance of Na background (mS/uF)

	// NSR Ca Ion Concentration Changes 
		var.kmup = 0.00092;   // Half-saturation concentration of iup (mM)
		var.iupbar = 0.00875; // Max. current through iup channel (mM/ms)
		var.nsrbar = 15;      // Max. [Ca] in NSR (mM)

	// JSR Ca Ion Concentration Changes 
		var.tauon = 2.0;        // Time constant of activation of Ca release from JSR (ms)
		var.tauoff = 2.0;       // Time constant of deactivation of Ca release from JSR (ms)
		var.csqnth = 8.75;    // Threshold for release of Ca from CSQN due to JSR ovreload (mM)
		var.gmaxrel = 150;    // Max. rate constant of Ca release from JSR due to overload (ms^-1)
		var.grelbarjsrol = 0; // Rate constant of Ca release from JSR due to overload (ms^-1)
		var.swspontan = 0;    // switch of spontaneous release
		var.csqnbar = 10;     // Max. [Ca] buffered in CSQN (mM)
		var.kmcsqn = 0.8;     // Equalibrium constant of buffering for CSQN (mM)

		var.csqn = 6.97978;
		var.trpn = 0.0143923;
		var.cmdn = 0.00257849;
	
		// Another parameter initial setting
		var.t_cicr = 25;
		var.t_overload = 25;
		var.boolien = 1;
		var.dICa_total = 0;

	// Translocation of Ca Ions from NSR to JSR
		var.tautr = 180;      // Time constant of Ca transfer from NSR to JSR (ms)

	// Myoplasmic Ca Ion Concentration Changes 
		var.cmdnbar = 0.050;   // Max. [Ca] buffered in CMDN (mM)
		var.trpnbar = 0.070;   // Max. [Ca] buffered in TRPN (mM)
		var.kmcmdn = 0.00238;  // Equalibrium constant of buffering for CMDN (mM)
		var.kmtrpn = 0.0005;   // Equalibrium constant of buffering for TRPN (mM)

	// Extracellular Ion Concentration Changes 
		var.tau_diff = 1000; // Diffusion Constant for Ion Movement from Bulk Medium to Cleft Space

}

