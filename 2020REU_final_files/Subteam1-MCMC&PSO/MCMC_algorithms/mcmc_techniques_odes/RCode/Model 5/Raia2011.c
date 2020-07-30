/* file Raia2011.c */
#include <R.h>
static double parms[17];
static double forc[2];

#define CD274mRNA_production parms[0]
#define DecoyR_binding parms[1]
#define JAK2_p_inhibition parms[2]
#define JAK2_phosphorylation parms[3]
#define Kon_IL13Rec parms[4]
#define Rec_intern parms[5]
#define Rec_phosphorylation parms[6]
#define Rec_recycle parms[7]
#define SOCS3_accumulation parms[8]
#define SOCS3_degradation parms[9]
#define SOCS3_translation parms[10]
#define SOCS3mRNA_production parms[11]
#define STAT5_phosphorylation parms[12]
#define pJAK2_dephosphorylation parms[13]
#define pRec_degradation parms[14]
#define pRec_intern parms[15]
#define pSTAT5_dephosphorylation parms[16]

/*

#define init_Rec_i parms[17] esta se debe calcular y derivar no encuentro este parámetro


Estos parámetros se utiliza en los observables para 
minimizar la función objeto
#define scale_CD274mRNA_obs parms[18]
#define scale_IL13_cell_obs parms[19]
#define scale_SOCS3mRNA_obs parms[20]
#define scale_pIL4Ra_obs parms[21]
#define scale_pJAK2_obs parms[22]
*/
#define SHP1 forc[0]
#define il13_level forc[1]

/*-------------------------*/

#define dRec ydot[0]
#define dRec_i ydot[1]
#define dIL13_Rec ydot[2]
#define dp_IL13_Rec ydot[3]
#define dp_IL13_Rec_i ydot[4]
#define dJAK2 ydot[5]
#define dpJAK2 ydot[6]
#define dSTAT5 ydot[7]
#define dpSTAT5 ydot[8]
#define dSOCS3mRNA ydot[9]
#define dDecoyR ydot[10]
#define dIL13_DecoyR ydot[11]
#define dSOCS3 ydot[12]
#define dCD274mRNA ydot[13]


/*xstart*/
#define Rec y[0]
#define Rec_i y[1]
#define IL13_Rec y[2]
#define p_IL13_Rec y[3]
#define p_IL13_Rec_i y[4]
#define JAK2 y[5]
#define pJAK2 y[6]
#define STAT5 y[7]
#define pSTAT5 y[8]
#define SOCS3mRNA y[9]
#define DecoyR y[10]
#define IL13_DecoyR y[11]
#define SOCS3 y[12]
#define CD274mRNA y[13]


/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=17;
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

	dRec = Rec_i*Rec_recycle - Rec*Rec_intern - il13_level*Kon_IL13Rec*Rec ;
	dRec_i = Rec*Rec_intern - Rec_i*Rec_recycle ;
	dIL13_Rec = il13_level*Kon_IL13Rec*Rec - IL13_Rec*Rec_phosphorylation*pJAK2 ;
	dp_IL13_Rec = IL13_Rec*Rec_phosphorylation*pJAK2 - pRec_intern*p_IL13_Rec ;
	dp_IL13_Rec_i = pRec_intern*p_IL13_Rec - pRec_degradation*p_IL13_Rec_i ;
	dJAK2 = SHP1*pJAK2*pJAK2_dephosphorylation - (JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1) ;
	dpJAK2 = (JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - SHP1*pJAK2*pJAK2_dephosphorylation + (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1) ;
	dSTAT5 = SHP1*pSTAT5*pSTAT5_dephosphorylation - STAT5*STAT5_phosphorylation*pJAK2 ;
	dpSTAT5 = STAT5*STAT5_phosphorylation*pJAK2 - SHP1*pSTAT5*pSTAT5_dephosphorylation ;
	dSOCS3mRNA = SOCS3mRNA_production*pSTAT5 ;
	dDecoyR = -il13_level*DecoyR*DecoyR_binding ;
	dIL13_DecoyR = il13_level*DecoyR*DecoyR_binding ;
	dSOCS3 = (SOCS3mRNA*SOCS3_translation)/(SOCS3mRNA + SOCS3_accumulation) - SOCS3*SOCS3_degradation ;
	dCD274mRNA = CD274mRNA_production*pSTAT5 ;

}

/*
c(45.3, 9.06, 181.2, 0)
Rec_i*Rec_recycle - Rec*Rec_intern - 2.265*Kon_IL13Rec*Rec*il13_level ;
Rec_i*Rec_recycle - Rec*Rec_intern - 9.06*Kon_IL13Rec*Rec
Rec_i*Rec_recycle - Rec*Rec_intern - 181.2*Kon_IL13Rec*Rec
Rec_i*Rec_recycle - Rec*Rec_intern
Rec_i*Rec_recycle - Rec*Rec_intern - 45.3*Kon_IL13Rec*Rec

Rec*Rec_intern - Rec_i*Rec_recycle 
Rec*Rec_intern - Rec_i*Rec_recycle
Rec*Rec_intern - Rec_i*Rec_recycle
Rec*Rec_intern - Rec_i*Rec_recycle
Rec*Rec_intern - Rec_i*Rec_recycle

2.265*Kon_IL13Rec*Rec*il13_level - IL13_Rec*Rec_phosphorylation*pJAK2
(3) 9.06*Kon_IL13Rec*Rec             - IL13_Rec*Rec_phosphorylation*pJAK2
(4) 181.2*Kon_IL13Rec*Rec            - IL13_Rec*Rec_phosphorylation*pJAK2
(1)                                  -IL13_Rec*Rec_phosphorylation*pJAK2
(2) 45.3*Kon_IL13Rec*Rec             - IL13_Rec*Rec_phosphorylation*pJAK2

IL13_Rec*Rec_phosphorylation*pJAK2 - pRec_intern*p_IL13_Rec 
IL13_Rec*Rec_phosphorylation*pJAK2 - pRec_intern*p_IL13_Rec
IL13_Rec*Rec_phosphorylation*pJAK2 - pRec_intern*p_IL13_Rec
IL13_Rec*Rec_phosphorylation*pJAK2 - pRec_intern*p_IL13_Rec
IL13_Rec*Rec_phosphorylation*pJAK2 - pRec_intern*p_IL13_Rec

pRec_intern*p_IL13_Rec - pRec_degradation*p_IL13_Rec_i 
pRec_intern*p_IL13_Rec - pRec_degradation*p_IL13_Rec_i
pRec_intern*p_IL13_Rec - pRec_degradation*p_IL13_Rec_i
pRec_intern*p_IL13_Rec - pRec_degradation*p_IL13_Rec_i
pRec_intern*p_IL13_Rec - pRec_degradation*p_IL13_Rec_i

SHP1*pJAK2*pJAK2_dephosphorylation - (JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1) ;
  91*pJAK2*pJAK2_dephosphorylation - (JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1)
  91*pJAK2*pJAK2_dephosphorylation - (JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1)
  91*pJAK2*pJAK2_dephosphorylation - (JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1)
  91*pJAK2*pJAK2_dephosphorylation - (JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1)

(JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - SHP1*pJAK2*pJAK2_dephosphorylation + (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1) ;
(JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - 91*pJAK2*pJAK2_dephosphorylation + (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1)
(JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - 91*pJAK2*pJAK2_dephosphorylation + (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1)
(JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - 91*pJAK2*pJAK2_dephosphorylation + (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1)
(JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - 91*pJAK2*pJAK2_dephosphorylation + (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1)

SHP1*pSTAT5*pSTAT5_dephosphorylation - STAT5*STAT5_phosphorylation*pJAK2 ;
91*pSTAT5*pSTAT5_dephosphorylation - STAT5*STAT5_phosphorylation*pJAK2
91*pSTAT5*pSTAT5_dephosphorylation - STAT5*STAT5_phosphorylation*pJAK2
91*pSTAT5*pSTAT5_dephosphorylation - STAT5*STAT5_phosphorylation*pJAK2
91*pSTAT5*pSTAT5_dephosphorylation - STAT5*STAT5_phosphorylation*pJAK2

STAT5*STAT5_phosphorylation*pJAK2 - SHP1*pSTAT5*pSTAT5_dephosphorylation ;
STAT5*STAT5_phosphorylation*pJAK2 - 91*pSTAT5*pSTAT5_dephosphorylation
STAT5*STAT5_phosphorylation*pJAK2 - 91*pSTAT5*pSTAT5_dephosphorylation
STAT5*STAT5_phosphorylation*pJAK2 - 91*pSTAT5*pSTAT5_dephosphorylation
STAT5*STAT5_phosphorylation*pJAK2 - 91*pSTAT5*pSTAT5_dephosphorylation

SOCS3mRNA_production*pSTAT5 ;
SOCS3mRNA_production*pSTAT5
SOCS3mRNA_production*pSTAT5
SOCS3mRNA_production*pSTAT5
SOCS3mRNA_production*pSTAT5

-2.265*DecoyR*DecoyR_binding*il13_level ;
-9.06*DecoyR*DecoyR_binding
-181.2*DecoyR*DecoyR_binding
0
-45.3*DecoyR*DecoyR_binding

2.265*DecoyR*DecoyR_binding*il13_level ;
9.06*DecoyR*DecoyR_binding
181.2*DecoyR*DecoyR_binding
0
45.3*DecoyR*DecoyR_binding

(SOCS3mRNA*SOCS3_translation)/(SOCS3mRNA + SOCS3_accumulation) - SOCS3*SOCS3_degradation 
(SOCS3mRNA*SOCS3_translation)/(SOCS3mRNA + SOCS3_accumulation) - SOCS3*SOCS3_degradation
(SOCS3mRNA*SOCS3_translation)/(SOCS3mRNA + SOCS3_accumulation) - SOCS3*SOCS3_degradation
(SOCS3mRNA*SOCS3_translation)/(SOCS3mRNA + SOCS3_accumulation) - SOCS3*SOCS3_degradation
(SOCS3mRNA*SOCS3_translation)/(SOCS3mRNA + SOCS3_accumulation) - SOCS3*SOCS3_degradation

CD274mRNA_production*pSTAT5
CD274mRNA_production*pSTAT5
CD274mRNA_production*pSTAT5
CD274mRNA_production*pSTAT5 
CD274mRNA_production*pSTAT5 */

