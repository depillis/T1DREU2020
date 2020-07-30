/* file smad.c */
#include <R.h>
static double parms[13];
#define KCAT parms[0]
#define K1 parms[1]
#define k5nc parms[2]
#define k5cn parms[3]
#define k4nc parms[4]
#define k4cn parms[5]
#define k2a parms[6]
#define k2d parms[7]
#define k3 parms[8]
#define k6d parms[9]
#define k6a parms[10]
#define Vmax7 parms[11]
#define K7 parms[12]

#define dreceptor ydot[0]
#define dR_smad_cyt ydot[1]
#define dR_smad_P_cyt ydot[2]
#define dsmad4_cyt ydot[3]
#define dR_smad_P_smad4_cyt ydot[4]
#define dR_smad_P_smad4_nuc ydot[5]
#define dR_smad_nuc ydot[6]
#define dR_smad_P_nuc ydot[7]
#define dsmad4_nuc ydot[8]
#define dPi ydot[9]

/*xstart*/
#define receptor y[0]
#define R_smad_cyt y[1]
#define R_smad_P_cyt y[2]
#define smad4_cyt y[3]
#define R_smad_P_smad4_cyt y[4]
#define R_smad_P_smad4_nuc y[5]
#define R_smad_nuc y[6]
#define R_smad_P_nuc y[7]
#define smad4_nuc y[8]
#define Pi y[9]


/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=13;
odeparms(&N, parms);
}
/* Derivatives and 1 output variable */


void derivs (int *neq, double *t, double *y, double *ydot,double *yout, int *ip)
{
	if (ip[0] <1) error("nout should be at least 1");  /* checking whether enough memory for output has been allocated*/


	/*dreceptor=-100*exp(-1/90*(*t));*/

	dreceptor=-1*receptor/90;

	/*dreceptor=-1*receptor/90;*/

	dR_smad_cyt=-KCAT*receptor*R_smad_cyt/(K1+R_smad_cyt)+k5nc*R_smad_nuc-k5cn*R_smad_cyt;

	dR_smad_P_cyt=KCAT*receptor*R_smad_cyt/(K1+R_smad_cyt)-k2a*R_smad_P_cyt*smad4_cyt+k2d*R_smad_P_smad4_cyt;

	dsmad4_cyt=k4nc*smad4_nuc-k4cn*smad4_cyt-k2a*R_smad_P_cyt*smad4_cyt+k2d*R_smad_P_smad4_cyt;

	dR_smad_P_smad4_cyt=k2a*R_smad_P_cyt*smad4_cyt-k2d*R_smad_P_smad4_cyt-k3*R_smad_P_smad4_cyt;

	dR_smad_P_smad4_nuc=k3*R_smad_P_smad4_cyt-k6d*R_smad_P_smad4_nuc+k6a*smad4_nuc*R_smad_P_nuc;

	dR_smad_nuc=-k5nc*R_smad_nuc+k5cn*R_smad_cyt+Vmax7*R_smad_P_nuc/(K7+R_smad_P_nuc);

	dR_smad_P_nuc=k6d*R_smad_P_smad4_nuc-k6a*smad4_nuc*R_smad_P_nuc-Vmax7*R_smad_P_nuc/(K7+R_smad_P_nuc);

	dsmad4_nuc=-k4nc*smad4_nuc+k4cn*smad4_cyt+k6d*R_smad_P_smad4_nuc-k6a*smad4_nuc*R_smad_P_nuc;

	dPi=Vmax7*R_smad_P_nuc/(K7+R_smad_P_nuc);


}

