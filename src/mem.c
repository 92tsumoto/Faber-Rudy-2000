#include "syspara.h"

typedef double Number;
typedef long long Lint;

void mem()
{
	int i,k,m,l,z;

	// initialized tablization memorys for Exp functions
	var.Tam=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbm=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tah=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbh=(Number *)calloc(VNMAX,sizeof(Number));
	var.Taj=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbj=(Number *)calloc(VNMAX,sizeof(Number));
	
	var.Tdss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttaud=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tfss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttauf=(Number *)calloc(VNMAX,sizeof(Number));
	
	var.Tbss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttaub=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tgss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttaug=(Number *)calloc(VNMAX,sizeof(Number));
	
	var.Txrss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttauxr=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tr=(Number *)calloc(VNMAX,sizeof(Number));
	var.Txr2ss=(Number *)calloc(VNMAX,sizeof(Number));
	
	var.Txs1ss=(Number *)calloc(VNMAX,sizeof(Number));
	var.Ttauxs1=(Number *)calloc(VNMAX,sizeof(Number));

	var.Tkp=(Number *)calloc(VNMAX,sizeof(Number));
	
	var.Tpov=(Number *)calloc(VNMAX,sizeof(Number));

	var.Trvdv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tazdv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbzdv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Taydv=(Number *)calloc(VNMAX,sizeof(Number));
	var.Tbydv=(Number *)calloc(VNMAX,sizeof(Number));

	var.TexpCa=(Number *)calloc(VNMAX,sizeof(Number));
	var.TexpNa=(Number *)calloc(VNMAX,sizeof(Number));
	var.TexpK=(Number *)calloc(VNMAX,sizeof(Number));

	var.Texp0=(Number *)calloc(VNMAX,sizeof(Number));
	var.Texp1=(Number *)calloc(VNMAX,sizeof(Number));
	var.Texp2=(Number *)calloc(VNMAX,sizeof(Number));
	var.Texp3=(Number *)calloc(VNMAX,sizeof(Number));
	if(var.Tam==NULL || var.Tah==NULL || var.Taj==NULL || var.Tbm==NULL || var.Tbh==NULL || var.Tbj==NULL
		|| var.Tdss==NULL || var.Ttaud==NULL || var.Tfss==NULL || var.Ttauf==NULL 
		|| var.TexpCa==NULL || var.TexpNa==NULL || var.TexpK==NULL 
		|| var.Txs1ss==NULL || var.Ttauxs1==NULL 
		|| var.Txrss==NULL || var.Ttauxr==NULL || var.Tr==NULL || var.Txr2ss==NULL
		|| var.Tkp==NULL 
		|| var.Tpov==NULL 
		|| var.Tbss==NULL || var.Ttaub==NULL || var.Tgss==NULL || var.Ttaug==NULL
		|| var.Trvdv==NULL || var.Tazdv==NULL || var.Tbzdv==NULL || var.Taydv==NULL || var.Tbydv==NULL 
		|| var.Texp0==NULL || var.Texp1==NULL 
		|| var.Texp2==NULL || var.Texp3==NULL 
	) exit(1);

}
		
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
void close_mem()
{

	free(var.Tam);free(var.Tah);free(var.Taj);free(var.Tbm);free(var.Tbh);free(var.Tbj);
	free(var.TexpCa);free(var.TEna);free(var.TEk);
	free(var.Tdss);free(var.Ttaud);free(var.Tfss);free(var.Ttauf);
	free(var.Tbss);free(var.Ttaub);free(var.Tgss);free(var.Ttaug);
	free(var.Txrss);free(var.Ttauxr);free(var.Tr);free(var.Txr2ss);
	free(var.Txs1ss);free(var.Ttauxs1);
	free(var.Tkp);
	free(var.Tpov);
	free(var.Trvdv);free(var.Tazdv);free(var.Tbzdv);free(var.Taydv);free(var.Tbydv);
	free(var.Texp0);free(var.Texp1);free(var.Texp2);free(var.Texp3);

}
