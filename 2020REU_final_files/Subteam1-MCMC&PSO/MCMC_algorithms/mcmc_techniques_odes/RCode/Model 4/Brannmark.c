/* file Brannmark.c */
#include <R.h>
static double parms[14];
static double forc[2];

#define k1a parms[0]
#define k1aBasic parms[1]
#define k1b parms[2]
#define k1c parms[3]
#define k1d parms[4]
#define k1e parms[5]
#define k1f parms[6]
#define k1g parms[7]
#define k1r parms[8]
#define k21 parms[9]
#define k22 parms[10]
#define k3 parms[11]
#define km2 parms[12]
#define km3 parms[13]
/*
Estos parámetros se utiliza en los observables para 
minimizar la función objeto
#define k_IRP_1Step parms[14]
#define k_IRSiP_1Step parms[15]
#define k_IRSiP_2Step parms[16]
#define k_IRSiP_DosR parms[17]
#define sigmaY1TimR parms[18]
#define sigmaY2Step parms[19]
#define sigmaY2TimR parms[20]
#define sigmaYDosR parms[21]
*/
#define insulin1 forc[0]
#define insulin2 forc[1]

/*-------------------------*/

#define dIR ydot[0]
#define dIRins ydot[1]
#define dIRp ydot[2]
#define dIRiP ydot[3]
#define dIRi ydot[4]
#define dIRS ydot[5]
#define dIRSiP ydot[6]
#define dXX ydot[7]
#define dXXp ydot[8]

/*xstart*/
#define IR y[0]
#define IRins y[1]
#define IRp y[2]
#define IRiP y[3]
#define IRi y[4]
#define IRS y[5]
#define IRSiP y[6]
#define XX y[7]
#define XXp y[8]



/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=14;
odeparms(&N, parms);
}
/* Derivatives and 1 output variable */


void forcc(void (* odeforcs)(int *, double *))
{
int N=2;
odeforcs(&N, forc);
}


void derivs (int *neq, double *t, double *y, double *ydot,double *yout, int *ip)
{
	if (ip[0] <1) error("nout should be at least 1");  /* checking whether enough memory for output has been allocated*/

	dIR = IRins*k1b - IR*k1aBasic + IRp*k1g + IRi*k1r - IR*(insulin2 + insulin1)*k1a; 
	dIRins = IR*k1aBasic - IRins*k1b - IRins*k1c + IR*(insulin2 + insulin1)*k1a;
	dIRp = IRins*k1c - IRp*k1d - IRp*k1g;
	dIRiP = IRp*k1d - IRiP*(k1e + (XXp*k1f)/(XXp + 1));
	dIRi = IRiP*(k1e + (XXp*k1f)/(XXp + 1)) - IRi*k1r;
	dIRS = IRSiP*km2 - IRS*k21*(IRp + IRiP*k22);
	dIRSiP = IRS*k21*(IRp + IRiP*k22) - IRSiP*km2;
	dXX = XXp*km3 - IRSiP*XX*k3;
	dXXp = IRSiP*XX*k3 - XXp*km3;
}
