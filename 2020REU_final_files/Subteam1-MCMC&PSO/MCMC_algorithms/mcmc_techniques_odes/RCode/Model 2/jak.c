/* file jak.c */
#include <R.h>
static double parms[26];
#define k1 parms[0]    	/*  JAK2_phosphorylation_by_Epo   		*/
#define k2 parms[1]    	/*  EpoR_phosphorylation_by_pJAK2 		*/
#define k3 parms[2]    	/*  SHP1_activation_by_pEpoR      		*/
#define k4 parms[3]    	/*  SHP1_delay			 		*/
#define k13 parms[4]  	/*  actSHP1_deactivation			*/
#define k14 parms[5]  	/*  pEpoR_desphosphorylation_by_actSHP1 	*/
#define k15 parms[6] 	/*  pJAK2_dephosphorylation_by_actSHP		*/
#define k16 parms[7]  	/*  SOS_recruitment_by_pEpoR			*/
#define k17 parms[8]  	/*  mSOS_release_from_membrane			*/
#define k18 parms[9]  	/*  mSOS_induced_Raf_phosphorylation		*/
#define k19 parms[10] 	/*  pRaf_dephosphorylation			*/
#define k20 parms[11]	/*  First_MEK2_phosphorylation_by_pRaf		*/
#define k21 parms[12]	/*  First_MEK1_phosphorylation_by_pRaf		*/
#define k22 parms[13]	/*  Second_MEK2_phosphorylation_by_pRaf		*/
#define k23 parms[14]	/*  Second_MEK1_phosphorylation_by_pRaf		*/
#define k24 parms[15]	/*  First_MEK_desphosphorylation		*/
#define k26 parms[16]	/*  Second_MEK_desphosphorylation		*/
#define k28 parms[17]	/*  First_ERK1_phosphorylation_by_ppMEK		*/
#define k29 parms[18]	/*  First_ERK2_phosphorylation_by_ppMEK		*/
#define k32 parms[19]	/*  Second_ERK1_phosphorylation_by_ppMEK	*/
#define k33 parms[20]	/*  Second_ERK2_phosphorylation_by_ppMEK	*/
#define k36 parms[21]	/*  First_ERK_dephosphorylation			*/
#define k38 parms[22]	/*  Second_ERK_dephosphorylation		*/
#define k40 parms[23]	/*  ppERK_neg_feedback_on_mSOS			*/
#define k42 parms[24]	/*  pSOS_dephosphorylation			*/
#define cell parms[25]

#define dJAK2 ydot[0]
#define dpJAK2 ydot[1]
#define dEpo ydot[2]
#define dEpoR ydot[3]
#define dpEpoR ydot[4]
#define dSHP1 ydot[5]
#define dmSHP1 ydot[6]
#define dDelay01_mSHP1 ydot[7]
#define dDelay02_mSHP1 ydot[8]
#define dDelay03_mSHP1 ydot[9]
#define dDelay04_mSHP1 ydot[10]
#define dDelay05_mSHP1 ydot[11]
#define dDelay06_mSHP1 ydot[12]
#define dDelay07_mSHP1 ydot[13]
#define dDelay08_mSHP1 ydot[14]
#define dactSHP1 ydot[15]
#define dSOS ydot[16]
#define dmSOS ydot[17]
#define dRaf ydot[18]
#define dpRaf ydot[19]
#define dMEK2 ydot[20]
#define dpMEK2 ydot[21]
#define dMEK1 ydot[22]
#define dpMEK1 ydot[23]
#define dppMEK2 ydot[24]
#define dppMEK1 ydot[25]
#define dERK1 ydot[26]
#define dpERK1 ydot[27]
#define dERK2 ydot[28]
#define dpERK2 ydot[29]
#define dppERK1 ydot[30]
#define dppERK2 ydot[31]
#define dpSOS ydot[32]

/*xstart*/
#define JAK2 y[0]
#define pJAK2 y[1]
#define Epo y[2]
#define EpoR y[3]
#define pEpoR y[4]
#define SHP1 y[5]
#define mSHP1 y[6]
#define Delay01_mSHP1 y[7]
#define Delay02_mSHP1 y[8]
#define Delay03_mSHP1 y[9]
#define Delay04_mSHP1 y[10]
#define Delay05_mSHP1 y[11]
#define Delay06_mSHP1 y[12]
#define Delay07_mSHP1 y[13]
#define Delay08_mSHP1 y[14]
#define actSHP1 y[15]
#define SOS y[16]
#define mSOS y[17]
#define Raf y[18]
#define pRaf y[19]
#define MEK2 y[20]
#define pMEK2 y[21]
#define MEK1 y[22]
#define pMEK1 y[23]
#define ppMEK2 y[24]
#define ppMEK1 y[25]
#define ERK1 y[26]
#define pERK1 y[27]
#define ERK2 y[28]
#define pERK2 y[29]
#define ppERK1 y[30]
#define ppERK2 y[31]
#define pSOS y[32]


/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=26;
odeparms(&N, parms);
}
/* Derivatives and 1 output variable */


void derivs (int *neq, double *t, double *y, double *ydot,double *yout, int *ip)
{
	if (ip[0] <1) error("nout should be at least 1");  /* checking whether enough memory for output has been allocated*/




	dJAK2= -k1*JAK2*Epo*cell+k15*pJAK2*actSHP1*cell;

	dpJAK2= k1*JAK2*Epo*cell-k15*pJAK2*actSHP1*cell;

	dEpo= 0;

	dEpoR= -k2*EpoR*pJAK2*cell+k14*pEpoR*actSHP1*cell;

	dpEpoR= k2*EpoR*pJAK2*cell-k14*pEpoR*actSHP1*cell;
	
	dSHP1= -k3*SHP1*pEpoR*cell+k13*actSHP1*cell;

	dmSHP1= k3*SHP1*pEpoR*cell-k4*mSHP1*cell;

	dDelay01_mSHP1= k4*mSHP1*cell-k4*Delay01_mSHP1*cell;

	dDelay02_mSHP1= k4*Delay01_mSHP1*cell-k4*Delay02_mSHP1*cell;

	dDelay03_mSHP1= k4*Delay02_mSHP1*cell-k4*Delay03_mSHP1*cell;

	dDelay04_mSHP1= k4*Delay03_mSHP1*cell-k4*Delay04_mSHP1*cell;

	dDelay05_mSHP1= k4*Delay04_mSHP1*cell-k4*Delay05_mSHP1*cell;

	dDelay06_mSHP1= k4*Delay05_mSHP1*cell-k4*Delay06_mSHP1*cell;	

	dDelay07_mSHP1= k4*Delay06_mSHP1*cell-k4*Delay07_mSHP1*cell;

	dDelay08_mSHP1= k4*Delay07_mSHP1*cell-k4*Delay08_mSHP1*cell;

	dactSHP1= k4*Delay08_mSHP1*cell-k13*actSHP1*cell;

	dSOS= -k16*SOS*pEpoR*cell+k17*mSOS*cell+k42*pSOS*cell;

	dmSOS= k16*SOS*pEpoR*cell-k17*mSOS*cell-k40*mSOS*ppERK1*cell-k40*mSOS*ppERK2*cell;
	
	dRaf= -k18*Raf*mSOS*cell+k19*pRaf*cell;
	
	dpRaf= k18*Raf*mSOS*cell-k19*pRaf*cell;

	dMEK2= -k20*MEK2*pRaf*cell+k26*pMEK2*cell;

	dpMEK2= k20*MEK2*pRaf*cell-k22*pMEK2*pRaf*cell+k24*ppMEK2*cell-k26*pMEK2*cell;

	dMEK1= -k21*MEK1*pRaf*cell+k26*pMEK1*cell;

	dpMEK1= k21*MEK1*pRaf*cell-k23*pMEK1*pRaf*cell+k24*ppMEK1*cell-k26*pMEK1*cell;

	dppMEK2= k22*pMEK2*pRaf*cell-k24*ppMEK2*cell;

	dppMEK1= k23*pMEK1*pRaf*cell-k24*ppMEK1*cell;

	dERK1= -k28*ERK1*ppMEK2*cell-k28*ERK1*ppMEK1*cell+k38*pERK1*cell;

	dpERK1= k28*ERK1*ppMEK2*cell+k28*ERK1*ppMEK1*cell-k32*pERK1*ppMEK2*cell-k32*pERK1*ppMEK1*cell+k36*ppERK1*cell-k38*pERK1*cell;

	dERK2= -k29*ERK2*ppMEK2*cell-k29*ERK2*ppMEK1*cell+k38*pERK2*cell;

	dpERK2= k29*ERK2*ppMEK2*cell+k29*ERK2*ppMEK1*cell-k33*pERK2*ppMEK2*cell-k33*pERK2*ppMEK1*cell+k36*ppERK2*cell-k38*pERK2*cell;

	dppERK1= k32*pERK1*ppMEK2*cell+k32*pERK1*ppMEK1*cell-k36*ppERK1*cell;

	dppERK2= k33*pERK2*ppMEK2*cell+k33*pERK2*ppMEK1*cell-k36*ppERK2*cell;

	dpSOS= k40*mSOS*ppERK1*cell+k40*mSOS*ppERK2*cell-k42*pSOS*cell;
	


}

