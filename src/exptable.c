#include "syspara.h"

void make_ExpTable()
{

	int vindex,naiindex,kiindex,tindex;
	double v,nai,ki,ton;
    
	for(vindex=0;vindex<VNMAX;vindex++){

		v = (double)vindex/dvm-(double)Emax;

		
        /** for ina **/
        var.Tam[vindex] = 0.32*(v+47.13)/(1.0-exp(-0.1*(v+47.13)));
        var.Tbm[vindex] = 0.08*exp(-v/11.0);
        
		if ( v >= -40.0 ) { // V>=40mV
            var.Tah[vindex] = 0.0;
            var.Tbh[vindex] = 1.0/(0.13*(1.0+exp((v+10.66)/-11.1)));
            var.Taj[vindex] = 0.0;
            var.Tbj[vindex] = (0.3*exp(-0.0000002535*v))/(1.0+exp(-0.1*(v+32.0)));
        } else { // V < -40mV
            var.Tah[vindex] = 0.135 * exp(-(80.0+v)/6.8);
            var.Tbh[vindex] = 3.56*exp(0.079*v)+3.1e5 * exp(0.35*v);
            var.Taj[vindex] = (-127140.0*exp(0.2444*v)-0.00003474*exp(-0.04391*v))*(v+37.78)/(1.0+exp(0.311*(v+79.23)));
            var.Tbj[vindex] = (0.1212*exp(-0.01052*v))/(1.0+exp(-0.1378*(v+40.14)));
        }

		// for ical
        var.Tdss[vindex] = 1.0/(1.0+exp(-(v+10.0)/6.24));
        if( fabs(v+10.0)<0.01 ) {
            var.Ttaud[vindex] = var.Tdss[vindex]*(25.0/156.0)/(0.035);
        } else {
            var.Ttaud[vindex] = var.Tdss[vindex]*(1.0-exp(-(v+10.0)/6.24))/(0.035*(v+10.0));
        }

        var.Tfss[vindex] = 1.0/(1.0+exp((v+32.0)/8.0))+(0.6/(1.0+exp((50.0-v)/20.0)));;
        var.Ttauf[vindex] = 1.0/(0.0197*exp(-pow(0.0337*(v+10.0),2))+0.02); 

        var.TexpCa[vindex] = exp(v/var.RTon2F);
        var.TexpNa[vindex] = exp(v/var.RTonF);
        var.TexpK[vindex] = exp(v/var.RTonF);
		
		// icat
        var.Tbss[vindex] = 1.0/(1.0+exp(-(v+14.0)/10.8));
        var.Ttaub[vindex] = 3.7+6.1/(1.0+exp((v+25.0)/4.5 ));
        var.Tgss[vindex] = 1.0/(1.0+exp((v+60.0)/5.6));
        if ( v <= 0 ) {
            var.Ttaug[vindex] = -0.875*v+12.0;
        } else {
            var.Ttaug[vindex] = 12.0;
        }

        // for ikr 
        var.Txrss[vindex] = 1.0/(1.0+exp(-(v+21.5)/7.5));
        var.Ttauxr[vindex] = 1.0/(0.00138*(v+14.2)/(1.0-exp(-0.123*(v+14.2)))+0.00061*(v+38.9)/(exp(0.145*(v+38.9))-1.0));
        var.Tr[vindex] = 1.0/(1.0+exp((v+9.0)/22.4));
        
		var.Txr2ss[vindex] = 1.0/(1.0+exp(-(v+21.5+26.5)/7.5));
        

        // for iks 

        var.Txs1ss[vindex] = 1.0/(1.0+exp(-(v-1.5)/16.7));

		if(fabs(v+30.0) > 0.00001) {    // NaN will be occured at v=-30
			var.Ttauxs1[vindex] = 1.0/(0.0000719*(v+30.0)/(1.0-exp(-0.148*(v+30.0)))+0.000131*(v+30.0)/(exp(0.0687*(v+30.0))-1.0));
	    } else {
			var.Ttauxs1[vindex] = 1.0 / (0.0000719/(0.148-0.5*0.148*0.148*(v+30.0))+0.000131/(0.0687+0.5*0.0687*0.0687*(v+30.0)));
		}

		// ikp 
        var.Tkp[vindex] = 1.0/(1.0+exp((7.488-v)/5.98));
	
		// ikna
		var.Tpov[vindex] = 0.8-0.65/(1.0+exp((v+125.0)/15.0));

		// ito
		var.Trvdv[vindex] = exp(v/100.0);
		var.Tazdv[vindex] = (10.0*exp((v-40.0)/25.0))/(1.0 + exp((v-40.0)/25.0));
		var.Tbzdv[vindex] = (10.0*exp(-(v+90.0)/25.0))/(1.0 + exp(-(v+90.0)/25.0));
		var.Taydv[vindex] = 0.015/(1.0+exp((v+60.0)/5.0));
		var.Tbydv[vindex] = (0.1*exp((v+25.0)/5.0))/(1.0+exp((v+25.0)/5.0));
        
		// inaca
        var.Texp2[vindex] = exp((var.gammas-1.0)*v/var.RTonF);
        var.Texp3[vindex] = exp(v/var.RTonF);

		// inak 
        var.Texp0[vindex] = 1.0+0.1245*exp(-0.1*v/var.RTonF);
        var.Texp1[vindex] = exp(-v/var.RTonF);

	}
	
	// Ena and Nak ATPase
	//for(naiindex=0;naiindex<NAIMAX;naiindex++){
	//	nai = (double)naiindex/dvm;
	//	var.TEna[naiindex] = var.RTonF*log(var.nao/nai);
	//	var.Tnak[naiindex] = 1.0/(1.0+pow(var.kmnai/nai,2.0));
	//}       
	// Ek
	//for(kiindex=0;kiindex<KIMAX;kiindex++){
	//	ki = (double)kiindex/dvm;
	//	var.TEk[kiindex] = var.RTonF*log(var.ko/ki);
		// for Ik1      
	//		var.T5Ki[kiindex]=exp(-0.2385*var.RTonF*log(var.ko/ki));
	//		var.T6Ki[kiindex]=exp(-0.08032*var.RTonF*log(var.ko/ki));
	//		var.T7Ki[kiindex]=exp(-0.06175*var.RTonF*log(var.ko/ki));
	//		var.T8Ki[kiindex]=exp(0.5143*var.RTonF*log(var.ko/ki));
	//}       



}
