/* file mapkk_feedback.c */
#include <R.h>
static double parms[22];
#define V1 parms[0]
#define n parms[1]
#define Ki parms[2]
#define K1 parms[3]
#define V2 parms[4]
#define K2 parms[5]
#define k3 parms[6]
#define K3 parms[7]
#define k4 parms[8]
#define K4 parms[9]
#define V5 parms[10]
#define K5 parms[11]
#define V6 parms[12]
#define K6 parms[13]
#define k7 parms[14]
#define K7 parms[15]
#define k8 parms[16]
#define K8 parms[17]
#define V9 parms[18]
#define K9 parms[19]
#define V10 parms[20]
#define K10 parms[21]

#define dMKKK  ydot[0]
#define dMKKK_P  ydot[1]
#define dMKK  ydot[2]
#define dMKK_P ydot[3]
#define dMKK_PP  ydot[4]
#define dMAPK  ydot[5]
#define dMAPK_P  ydot[6]
#define dMAPK_PP  ydot[7]

/*xstart*/
#define MKKK  y[0]
#define MKKK_P  y[1]
#define MKK  y[2]
#define MKK_P y[3]
#define MKK_PP  y[4]
#define MAPK  y[5]
#define MAPK_P  y[6]
#define MAPK_PP  y[7]


/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=22;
odeparms(&N, parms);
}
/* Derivatives and 1 output variable */

void derivs (int *neq, double *t, double *y, double *ydot,double *yout, int *ip)
{
	if (ip[0] <1) error("nout should be at least 1");  /* checking whether enough memory for output has been allocated*/


 
	dMKKK= V2*MKKK_P/(K2+MKKK_P)-V1*MKKK/((1+pow((MAPK_PP/Ki),n))*(K1+MKKK));
	dMKKK_P= V1*MKKK/((1+pow((MAPK_PP/Ki),n))*(K1+MKKK)) - V2*MKKK_P/(K2+MKKK_P);	


	dMKK= V6*MKK_P/(K6+MKK_P) - k3*MKKK_P*MKK/(K3+MKK);
	dMKK_P= k3*MKKK_P*MKK/(K3+MKK) + V5*MKK_PP/(K5+MKK_PP) - k4*MKKK_P*MKK_P/(K4+MKK_P)  - V6*MKK_P/(K6+MKK_P);

	dMKK_PP= k4*MKKK_P*MKK_P/(K4+MKK_P)  - V5*MKK_PP/(K5+MKK_PP);
	dMAPK= V10*MAPK_P/(K10+MAPK_P)  - k7*MKK_PP*MAPK/(K7+MAPK);

	dMAPK_P= k7*MKK_PP*MAPK/(K7+MAPK) + V9*MAPK_PP/(K9+MAPK_PP) - k8*MKK_PP*MAPK_P/(K8+MAPK_P) - V10*MAPK_P/(K10+MAPK_P);
	dMAPK_PP= k8*MKK_PP*MAPK_P/(K8+MAPK_P)  - V9*MAPK_PP/(K9+MAPK_PP);

}

