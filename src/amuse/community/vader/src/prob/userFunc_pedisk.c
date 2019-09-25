#include "userFunc.h"
#include "math.h"

#include<stdio.h>

/**************************************************************************/
/* This defines userFunc routines for a protoplanetary disk that is being */
/* photoevaporated both by X-ray radiation from the central star and FUV  */
/* radiation by stars in the environment.                                 */
/* Parameter values (and their units) include:							  */
/* dMf_inner:	Internal photoevaporation mass loss rate (g/s)			  */
/* dMf_outer: 	External photoevaporation mass loss rate (g/s)			  */
/* C0: 			Column density at edge (g/cm^3)							  */
/* Tm: 	  		Midplane temperature (K)								  */
/* mmw:			Mean molecular weight (g)								  */
/* Mdot_nominal:Nominal accretion rate through inner boundary (g/s)		  */
/* 		WARNING: PAY ATTENTION TO THE UNITS, THIS IS NOT TRIVIAL		  */
/*																		  */
/* Written by Martijn Wilhelm, 2018										  */
/**************************************************************************/

void
userAlpha(const double t, const double dt, const grid *grd, 
	  const double *col, const double *pres, const double *eInt,
	  const double *gamma, const double *delta,
	  void *params,
	  double *alpha) {
  fprintf(stderr, 
	  "Warning: userAlpha function called but not implemented!\n");
  return;
}

void
userEOS(const double t, const double dt, const grid *grd, 
	const double *col, const double *pres, const double *eInt,
	void *params,
	double *gamma, double *delta) {
  fprintf(stderr, 
	  "Warning: userEOS function called but not implemented!\n");
  return;
}

void
userMassSrc(const double t, const double dt, const grid *grd,
	    const double *col, const double *pres, const double *eInt,
	    const double *gamma, const double *delta,
	    void *params,
	    double *massSrc) {
  fprintf(stderr, 
	  "Warning: userMassSrc function called but not implemented!\n");
  return;
}

void
userIntEnSrc(const double t, const double dt, const grid *grd,
	     const double *col, const double *pres, const double *eInt,
	     const double *gamma, const double *delta,
	     void *params, 
	     double *intEnSrc) {
  fprintf(stderr, 
	  "Warning: userIntEnSrc function called but not implemented!\n");
  return;
}

void
userIBC(const double t, const double dt, const grid *grd,
	const double *col, const double *pres, const double *eInt,
	const double *gamma, const double *delta,
	const pres_bc_type ibc_pres, const enth_bc_type ibc_enth,
	void *params, 
	double *ibc_pres_val, double *ibc_enth_val) {

  double C0			  = ((double *) params)[2];
  double Mdot_nominal = ((double *) params)[5];

  *ibc_enth_val = 0.;

  if (col[0]*grd->area[0] > Mdot_nominal*dt) {
	*ibc_pres_val = -Mdot_nominal;
  }
  else if (col[0] > C0) {
	*ibc_pres_val = -(col[0] - C0)/dt;
  }
  else { *ibc_pres_val = 0.; }

  return;
}

void
userOBC(const double t, const double dt, const grid *grd,
	const double *col, const double *pres, const double *eInt,
	const double *gamma, const double *delta,
	const pres_bc_type obc_pres, const enth_bc_type obc_enth,
	void *params, 
	double *obc_pres_val, double *obc_enth_val) {

  fprintf(stderr, 
	  "Warning: userOBC function called but not implemented!\n");
  return;
}

void
userPreTimestep(const double t, const double dt,
		const grid *grd, double *col, double *pres,
		double *eInt, double *mBnd, double *eBnd,
		double *mSrc, double *eSrc,
		void *params, const unsigned long nUserOut,
		double *userOut) {
  fprintf(stderr,
	  "Warning: userPreTimestep function called but not implemented!\n");

  return;
}

void
userPostTimestep(const double t, const double dt,
		 const grid *grd, double *col, double *pres,
		 double *eInt, double *mBnd, double *eBnd,
		 double *mSrc, double *eSrc,
		 void *params, const unsigned long nUserOut,
		 double *userOut) {

  double dMf_inner	= ((double *) params)[0];
  double dMf_outer	= ((double *) params)[1];
  double C0			= ((double *) params)[2];
  double Tm			= ((double *) params)[3];
  double mmw		= ((double *) params)[4];

  double cm_in_AU = 1.49597870691e+13; //according to amuse
  double kB = 1.3806504e-16; // erg/K, according to amuse
  double G = 6.674e-8; //cm3/g/s2

  double dM_inner = dMf_inner*dt;
  double dM_outer = dMf_outer*dt;
  double dM_inner_out = 0.0;
  double dM_outer_out = 0.0;
  double dM = 0.0;

  int edge = grd->nr-1;

  for (int i = grd->nr-1; i >= 0; i--) {
    if (col[i] > C0) { edge = i+1; break; }
  }
  if (edge == grd->nr-1) { edge -= 1; }


  //internal photoevaporation
  for (int i = 0; (i < grd->nr && dM_inner_out < dM_inner);
		i++) {
    dM = (col[i] - C0)*grd->area[i];

    if (dM > dM_inner - dM_inner_out) {
      col[i] -= (dM_inner - dM_inner_out)/grd->area[i];
      dM_inner_out = dM_inner;
    }
    else {
      col[i] -= dM/grd->area[i];
      dM_inner_out += dM;
    }
  }

  //external photoevaporation
  for (int i = edge; (i >= 0 && dM_outer_out < dM_outer); 
		i--) {
	dM = (col[i] - C0)*grd->area[i];

    if (dM > dM_outer - dM_outer_out) {
      col[i] -= (dM_outer - dM_outer_out)/grd->area[i];
      dM_outer_out = dM_outer;
    }
    else {
      col[i] -= dM/grd->area[i];
      dM_outer_out += dM;
    }
  }

  double T;

  //correct pressure
  for (int i = 0; i < grd->nr; i++) {
    T = Tm/sqrt(grd->r_g[i+1]/cm_in_AU);
    pres[i] = col[i]/mmw*kB*T;
  }

  return;
}

void
userCheckRead(
	      FILE *fp, grid *grd, const unsigned long nOut,
	      double *tOut, double *colOut,
	      double *presOut, double *eIntOut, double *mBndOut,
	      double *eBndOut, double *mSrcOut, double *eSrcOut,
	      const unsigned long nUserOut, double *userOut,
	      void *params
	      ) {
  fprintf(stderr,
	  "Warning: userCheckRead function called but not implemented!\n");
  return;
}

void
userCheckWrite(
	      FILE *fp,
	      const grid *grd, const unsigned long nOut,
	      const double *tOut, const double *colOut,
	      const double *presOut, const double *eIntOut,
	      const double *mBndOut, const double *eBndOut,
	      const double *mSrcOut, const double *eSrcOut,
	      const unsigned long nUserOut, const double *userOut,
	      const void *params
	      ) {
  fprintf(stderr,
	  "Warning: userCheckWrite function called but not implemented!\n");
  return;
}
