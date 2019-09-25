extern "C" {
	#include"src/vader_common.h"
	#include"src/driver.h"
	#include"src/init.h"
}

#include<iostream>
using namespace std;

/*
 * Interface code
 */

class outputs {
	public:

	 outputs();
	~outputs();

	unsigned long *nStep;
	unsigned long *nIter;
	unsigned long *nFail;

	unsigned long *nOut;
	double *tOut;
	double *colOut;
	double *presOut;
	double *eIntOut;
	double *mBndOut;
	double *eBndOut;
	double *mSrcOut;
	double *eSrcOut;
	double *userOut;
};

class control_params {
	public:

	double dtStart;
	double dtMin;
	double dtTol;
	double errTol;
	double maxDtIncrease;
	unsigned long maxIter;
	unsigned long interpOrder;
	unsigned long maxStep;
	bool useBE;
	bool preTimestep_func;
	bool postTimestep_func;
	unsigned long verbosity;
};

struct boundary_conditions {
	public:
	pres_bc_type	ibc_pres;
	double 			ibc_pres_val;
	enth_bc_type	ibc_enth;
	double			ibc_enth_val;
	bool			ibc_func;

	pres_bc_type	obc_pres;
	double 			obc_pres_val;
	enth_bc_type	obc_enth;
	double			obc_enth_val;
	bool			obc_func;
};

class Disk {
	public:

	 Disk();
	~Disk();

	grid *grd;
	wksp *w;
	boundary_conditions *bc;
	outputs *out;
	control_params *ctrl;

	int n;

	double current_time;

	double *column_density;
	double *pressure;
	double *internal_energy;

	bool alpha_func;
	double alpha;

	bool eos;
	double gamma;
	double delta;

	bool mass_source_func;
	double mass_source_val;

	bool internal_energy_source_func;
	double internal_energy_source_val;

	int number_of_parameters;
	double *params;
};

Disk::Disk() {

	bc = new boundary_conditions;
	ctrl = new control_params;	
	current_time = 0.;

	out = NULL;
	params = NULL;

	grd = NULL;
	w = NULL;

	column_density = NULL;
	pressure = NULL;
	internal_energy = NULL;
}

Disk::~Disk() {

	delete bc;

	delete column_density;
	delete pressure;
	delete internal_energy;

	delete params;
	delete out;
	delete ctrl;

	if (grd != NULL) { gridFree(grd); }
	else { delete grd; }	

	if (w != NULL) { wkspFree(w); }
	else { delete w; }

}

static Disk *disk;

outputs::outputs() {
	 nStep = new unsigned long;
	*nStep = 1;
	 nIter = new unsigned long;
	*nIter = 0;
	 nFail = new unsigned long;
	*nFail = 0;

	 nOut = new unsigned long;
	*nOut = 0;

	tOut    = new double[1];
	colOut  = new double[disk->n];
	presOut = new double[disk->n];
	eIntOut = new double[disk->n];
	mBndOut = new double[2];
		mBndOut[0] = 0.;
		mBndOut[1] = 0.;
	eBndOut = new double[2];
		eBndOut[0] = 0.;
		eBndOut[1] = 0.;
	mSrcOut = new double[disk->n];
	eSrcOut = new double[disk->n];
	userOut = new double[1];
}

outputs::~outputs() {

	delete nStep;
	delete nIter;
	delete nFail;
	delete nOut;
	delete [] tOut;
	delete [] colOut;
	delete [] presOut;
	delete [] eIntOut;
	delete [] mBndOut;
	delete [] eBndOut;
	delete [] mSrcOut;
	delete [] eSrcOut;
	delete [] userOut;

}


int initialize_code() { disk = new Disk; return 0; }


int initialize_keplerian_grid(const int n, const int linear, 
							  const double rmin, 
							  const double rmax, 
							  const double m) {

	delete [] disk->column_density;
	delete [] disk->pressure;
	delete [] disk->internal_energy;

	if (disk->grd != NULL) { gridFree(disk->grd); }
	else { delete disk->grd; }	

	if (disk->w != NULL) { wkspFree(disk->w); }
	else { delete disk->w; }

	delete disk->out;

	unsigned long N = static_cast<unsigned long>(n);
	int N_int = static_cast<int>(N);

	disk->grd = gridInitKeplerian(N, rmin, rmax, m, linear);
	disk->w = wkspAlloc(N);

	disk->column_density  = new double[N];
	disk->pressure 		  = new double[N];
	disk->internal_energy = new double[N];

	disk->n = N_int;
	disk->out = new outputs;

	return 0;
}

int update_keplerian_grid(const double m) {

	if (disk->grd == NULL) { return -1; }

	double rmin = disk->grd->r_h[0];
	double rmax = disk->grd->r_h[disk->n];
	unsigned int linear = disk->grd->linear;

	if (disk->grd != NULL) { gridFree(disk->grd); }
	else { delete disk->grd; }	

	if (disk->w != NULL) { wkspFree(disk->w); }
	else { delete disk->w; }

	disk->grd = gridInitKeplerian(disk->n, rmin, rmax, m, 
									linear);
	disk->w = wkspAlloc(disk->n);

	return 0;
}



int initialize_flat_grid(const int n, const int linear, 
						 const double rmin, const double rmax,
						 const double vphi) {

	delete [] disk->column_density;
	delete [] disk->pressure;
	delete [] disk->internal_energy;

	if (disk->grd != NULL) { gridFree(disk->grd); }
	else { delete disk->grd; }	

	if (disk->w != NULL) { wkspFree(disk->w); }
	else { delete disk->w; }

	delete disk->out;

	unsigned long N = static_cast<unsigned long>(n);
	int N_int = static_cast<int>(N);

	disk->grd = gridInitFlat(N, rmin, rmax, vphi, linear);
	disk->w = wkspAlloc(N);

	disk->column_density  = new double[N];
	disk->pressure 		  = new double[N];
	disk->internal_energy = new double[N];

	disk->n = N_int;
	disk->out = new outputs;

	return 0;
}

int update_flat_grid(const double vphi) {

	if (disk->grd == NULL) { return -1; }

	double rmin = disk->grd->r_h[0];
	double rmax = disk->grd->r_h[disk->n];
	unsigned int linear = disk->grd->linear;

	if (disk->grd != NULL) { gridFree(disk->grd); }
	else { delete disk->grd; }	

	if (disk->w != NULL) { wkspFree(disk->w); }
	else { delete disk->w; }

	disk->grd = gridInitFlat(disk->n, rmin, rmax, vphi, 
									linear);
	disk->w = wkspAlloc(disk->n);

	return 0;
}

/*
int initialize_tabulated_grid(const int n, const int linear, 
							  const double rmin, 
							  const double rmax, 
							  const int bspline_degree, 
							  const int bspline_breakpoints, 
							  const double *vTab, 
							  const double *rTab, 
							  const int nTab) {

	unsigned long N = static_cast<unsigned long>(n);

	disk->grd = gridInitTabulated(N, rTab, vTab, nTab, 
								  rmin, rmax, bspline_degree,
								  bspline_breakpoints, linear);
	disk->w = wkspAlloc(N);
	disk->n = N;

	disk->column_density  = new double[N];
	disk->pressure 		  = new double[N];
	disk->internal_energy = new double[N];

	return 0;
}

int initialize_grid(const int n, const int linear, 
					const double *r_g, const double *r_h,
					const double *vphi_g, const double *vphi_h, 
					const double *beta_g, const double *beta_h, 
					const double *psiEff_g, 
					const double *psiEff_h, 
					const double *g_h) {

	disk->grd = gridInitKeplerian(n, r_g, r_h, vphi_g, vphi_h,
								  beta_g, beta_h, 
								  psiEff_g, psiEff_h, 
								  g_h, linear);

	return 0;
}
*/


int evolve_model(double tlim) {

	double t_end;

	disk->out->tOut[0] = tlim;
	*(disk->out->nOut) = 0;

	bool *dummy_bool_pointer = new bool[1];
	dummy_bool_pointer[0] = false;

	t_end = driver(
			// Time parameters
		disk->current_time,
		tlim,
			// Equation of state parameters
		disk->eos,
		disk->gamma,
		disk->delta,
			// Dimensionless viscosity parameters
		disk->alpha_func,
		disk->alpha,
			// Inner boundary condition parameters
		disk->bc->ibc_pres,
		disk->bc->ibc_enth,
		disk->bc->ibc_func,
		disk->bc->ibc_pres_val,
		disk->bc->ibc_enth_val,
			// Outer boundary condition parameters
		disk->bc->obc_pres,
		disk->bc->obc_enth,
		disk->bc->obc_func,
		disk->bc->obc_pres_val,
		disk->bc->obc_enth_val,
			// Source function parameters
		disk->mass_source_func,
		disk->mass_source_val,
		disk->internal_energy_source_func,
		disk->internal_energy_source_val,
			// Control and method parameters
		disk->ctrl->dtStart,
		disk->ctrl->dtMin,
		disk->ctrl->dtTol,
		disk->ctrl->errTol,
		disk->ctrl->maxDtIncrease,
		disk->ctrl->maxIter,
		disk->ctrl->interpOrder,
		disk->ctrl->maxStep,
		disk->ctrl->useBE,
		disk->ctrl->preTimestep_func,
		disk->ctrl->postTimestep_func,
		disk->ctrl->verbosity,
			// Output control parameters
		1, disk->out->tOut, 0, 
		dummy_bool_pointer, dummy_bool_pointer, 
		"", false, false, 0,
			// Computational grid and workspace
		disk->grd,
		disk->w,
			// Input data
		disk->column_density,
		disk->pressure,
		disk->internal_energy,
			// User-defined extra parameters
		disk->params,
			// Diagnostics outputs
		disk->out->nStep,
		disk->out->nIter,
		disk->out->nFail,
			// Storage for Outputs
		disk->out->nOut,
		disk->out->tOut,
		disk->out->colOut,
		disk->out->presOut,
		disk->out->eIntOut,
		disk->out->mBndOut,
		disk->out->eBndOut,
		disk->out->mSrcOut,
		disk->out->eSrcOut,
		disk->out->userOut
	);

	disk->current_time = t_end;
	delete dummy_bool_pointer;

	if (t_end < tlim) { return -1; }
	else { return 0; }
}


int get_position_of_index(int i, double *r) {
	if (i < disk->n) {
		*r = disk->grd->r_g[i+1];
		return 0;
	}

	return -1;
}

int get_index_of_position(double r, int *i) {

	for (int j = 0; j < disk->n; j++) {
		if (disk->grd->r_g[j+1] == r) {
			*i = j;
			return 0;
		}
	}

	return -1;
}

int get_area_of_index(int *i, double *area, 
					  int number_of_points) {

	int i0;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			area[j] = disk->grd->area[i0];
		}
		else { return -1; }
	}

	return 0;
}

int get_rotational_velocity_of_index(int *i, double *vphi, 
									 int number_of_points) {
	
	int i0;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			vphi[j] = disk->grd->vphi_g[i0+1];
		}
		else { return -1; }
	}

	return 0;
}

int get_effective_potential_of_index(int *i, double *psiEff, 
									 int number_of_points) {
	
	int i0;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			psiEff[j] = disk->grd->psiEff_g[i0+1];
		}
		else { return -1; }
	}

	return 0;
}

int get_gravitational_potential_of_index(int *i, 
					double *psi_grav, int number_of_points) {
	
	int i0;
	double Ekin;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			Ekin = disk->grd->vphi_g[i0+1]*
				   disk->grd->vphi_g[i0+1]/2.;

			psi_grav[j] = disk->grd->psiEff_g[i0+1] - Ekin;
		}
		else { return -1; }
	}

	return 0;
}

int get_grid_column_density(int *i, double *column_density, 
							int number_of_points) {

	int i0;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			column_density[j] = disk->column_density[i0];
		}
		else { return -1; }
	}

	return 0;
}

int get_grid_pressure(int *i, double *pressure, 
					  int number_of_points) {

	int i0;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			pressure[j] = disk->pressure[i0];
		}
		else { return -1; }
	}

	return 0;
}

int get_grid_internal_energy(int *i, double *internal_energy, 
							 int number_of_points) {

	int i0;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			internal_energy[j] = disk->internal_energy[i0];
		}
		else { return -1; }
	}

	return 0;
}

int get_grid_state(int *i, double *column_density, 
				   double *pressure, double *internal_energy, 
				   int number_of_points) {

	int i0;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			column_density[j]  = disk->column_density[i0];
			pressure[j] 	   = disk->pressure[i0];
			internal_energy[j] = disk->internal_energy[i0];
		}
		else { return -1; }
	}

	return 0;
}

int set_grid_column_density(int *i, double *column_density, 
							int number_of_points) {

	int i0;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			disk->column_density[i0] = column_density[j];
		}
		else { return -1; }
	}

	return 0;
}

int set_grid_pressure(int *i, double *pressure, 
					  int number_of_points) {

	int i0;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			disk->pressure[i0] = pressure[j];
		}
		else { return -1; }
	}

	return 0;
}

int set_grid_internal_energy(int *i, double *internal_energy, 
							 int number_of_points) {

	int i0;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			disk->internal_energy[i0] = internal_energy[j];
		}
		else { return -1; }
	}

	return 0;
}

int set_grid_state(int *i, double *column_density, 
				   double *pressure, double *internal_energy, 
				   int number_of_points) {

	int i0;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			disk->column_density[i0]  = column_density[j];
			disk->pressure[i0] 		  = pressure[j];
			disk->internal_energy[i0] = internal_energy[j];
		}
		else { return -1; }
	}

	return 0;
}



int get_alpha_function(int *alpha_func) 
{ *alpha_func = disk->alpha_func; return 0; }

int set_alpha_function(int  alpha_func) 
{ disk->alpha_func = alpha_func;  return 0; }

int get_alpha(double *alpha) { *alpha = disk->alpha; return 0; }

int set_alpha(double  alpha) { disk->alpha = alpha;  return 0; }



int get_eos_function(int *eos) { *eos = disk->eos; return 0; }

int set_eos_function(int  eos) { disk->eos = eos;  return 0; }


int get_gamma(double *gamma) { *gamma = disk->gamma; return 0; }

int set_gamma(double  gamma) { disk->gamma = gamma;  return 0; }


int get_delta(double *delta) { *delta = disk->delta; return 0; }

int set_delta(double  delta) { disk->delta = delta;  return 0; }



int get_mass_source_function(int *mass_source_func) {
	*mass_source_func = disk->mass_source_func;  
	return 0; 
}

int set_mass_source_function(int  mass_source_func) { 
	disk->mass_source_func = mass_source_func;  
	return 0; 
}

int get_internal_energy_source_function
	(int *internal_energy_source_func) { 
	*internal_energy_source_func = 
		disk->internal_energy_source_func;
	return 0; 
}

int set_internal_energy_source_function
	(int  internal_energy_source_func) { 
	disk->internal_energy_source_func =
		internal_energy_source_func;  
	return 0; 
}


int get_mass_source_value(double *mass_source_val) 
{ *mass_source_val = disk->mass_source_val; return 0; }

int set_mass_source_value(double  mass_source_val)
{ disk->mass_source_val = mass_source_val;  return 0; }

int get_internal_energy_source_value
	(double *internal_energy_source_val) { 
	*internal_energy_source_val = 
		disk->internal_energy_source_val; 
	return 0; 
}

int set_internal_energy_source_value
	(double internal_energy_source_val) { 
	disk->internal_energy_source_val =
		internal_energy_source_val;  
	return 0; 
}



int get_inner_pressure_boundary_type(int *ibc_pres) {
	if	   (disk->bc->ibc_pres == FIXED_MASS_FLUX) 
	{ *ibc_pres = 1; return 0; }
	else if(disk->bc->ibc_pres == FIXED_TORQUE_FLUX) 
	{ *ibc_pres = 2; return 0; }
	else if(disk->bc->ibc_pres == FIXED_TORQUE)
	{ *ibc_pres = 3; return 0; }
	else { return -1; }
}

int set_inner_pressure_boundary_type(int  ibc_pres) {
	if	   (ibc_pres == 1) 
	{ disk->bc->ibc_pres = FIXED_MASS_FLUX;   return 0; }
	else if(ibc_pres == 2) 
	{ disk->bc->ibc_pres = FIXED_TORQUE_FLUX; return 0; }
	else if(ibc_pres == 3)
	{ disk->bc->ibc_pres = FIXED_TORQUE; 	  return 0; }
	else { return -1; }
}

int get_inner_pressure_boundary_mass_flux(double *ibc_pres_val)
{ *ibc_pres_val = disk->bc->ibc_pres_val; return 0; }

int set_inner_pressure_boundary_mass_flux(double  ibc_pres_val)
{  disk->bc->ibc_pres_val = ibc_pres_val; return 0; }

int get_inner_pressure_boundary_torque_flux(double *ibc_pres_val)
{ *ibc_pres_val = disk->bc->ibc_pres_val; return 0; }

int set_inner_pressure_boundary_torque_flux(double  ibc_pres_val)
{ disk->bc->ibc_pres_val = ibc_pres_val;  return 0; }

int get_inner_pressure_boundary_torque(double *ibc_pres_val)
{ *ibc_pres_val = disk->bc->ibc_pres_val; return 0; }

int set_inner_pressure_boundary_torque(double  ibc_pres_val)
{ disk->bc->ibc_pres_val = ibc_pres_val;  return 0; }


int get_inner_enthalpy_boundary_type(int *ibc_enth) { 
	if	   (disk->bc->ibc_enth == FIXED_ENTHALPY_VALUE) 
	{ *ibc_enth = 1; return 0; }
	else if(disk->bc->ibc_enth == FIXED_ENTHALPY_GRADIENT) 
	{ *ibc_enth = 2; return 0; }
	else { return -1; }
}

int set_inner_enthalpy_boundary_type(int ibc_enth) {
	if	   (ibc_enth == 1)
	{ disk->bc->ibc_enth = FIXED_ENTHALPY_VALUE;    return 0; }
	else if(ibc_enth == 2)
	{ disk->bc->ibc_enth = FIXED_ENTHALPY_GRADIENT; return 0; }
	else { return -1; }
}

int get_inner_enthalpy_boundary_enthalpy(double *ibc_enth_val)
{ *ibc_enth_val = disk->bc->ibc_enth_val; return 0; }

int set_inner_enthalpy_boundary_enthalpy(double  ibc_enth_val)
{ disk->bc->ibc_enth_val = ibc_enth_val;  return 0; }

int get_inner_enthalpy_boundary_enthalpy_gradient
	(double *ibc_enth_val)
{ *ibc_enth_val = disk->bc->ibc_enth_val; return 0; }

int set_inner_enthalpy_boundary_enthalpy_gradient
	(double  ibc_enth_val)
{ disk->bc->ibc_enth_val = ibc_enth_val;  return 0; }

int get_inner_boundary_function(int *ibc_func) 
{ *ibc_func = disk->bc->ibc_func; return 0; }

int set_inner_boundary_function(int  ibc_func) 
{ disk->bc->ibc_func = ibc_func;  return 0; }


int get_outer_pressure_boundary_type(int *obc_pres) {
	if	   (disk->bc->obc_pres == FIXED_MASS_FLUX) 
	{ *obc_pres = 1; return 0; }
	else if(disk->bc->obc_pres == FIXED_TORQUE_FLUX) 
	{ *obc_pres = 2; return 0; }
	else if(disk->bc->obc_pres == FIXED_TORQUE)
	{ *obc_pres = 3; return 0; }
	else { return -1; }
}

int set_outer_pressure_boundary_type(int obc_pres) {
	if	   (obc_pres == 1) 
	{ disk->bc->obc_pres = FIXED_MASS_FLUX;   return 0; }
	else if(obc_pres == 2) 
	{ disk->bc->obc_pres = FIXED_TORQUE_FLUX; return 0; }
	else if(obc_pres == 3)
	{ disk->bc->obc_pres = FIXED_TORQUE; 	  return 0; }
	else { return -1; }
}

int get_outer_pressure_boundary_mass_flux(double *obc_pres_val)
{ *obc_pres_val = disk->bc->obc_pres_val; return 0; }

int set_outer_pressure_boundary_mass_flux(double  obc_pres_val)
{ disk->bc->obc_pres_val = obc_pres_val;  return 0; }

int get_outer_pressure_boundary_torque_flux(double *obc_pres_val)
{ *obc_pres_val = disk->bc->obc_pres_val; return 0; }

int set_outer_pressure_boundary_torque_flux(double  obc_pres_val)
{ disk->bc->obc_pres_val = obc_pres_val;  return 0; }

int get_outer_pressure_boundary_torque(double *obc_pres_val)
{ *obc_pres_val = disk->bc->obc_pres_val; return 0; }

int set_outer_pressure_boundary_torque(double  obc_pres_val)
{ disk->bc->obc_pres_val = obc_pres_val;  return 0; }


int get_outer_enthalpy_boundary_type(int *obc_enth) { 
	if	   (disk->bc->obc_enth == FIXED_ENTHALPY_VALUE) 
	{ *obc_enth = 1; return 0; }
	else if(disk->bc->obc_enth == FIXED_ENTHALPY_GRADIENT) 
	{ *obc_enth = 2; return 0; }
	else { return -1; }
}

int set_outer_enthalpy_boundary_type(int obc_enth) {
	if	   (obc_enth == 1)
	{ disk->bc->obc_enth = FIXED_ENTHALPY_VALUE;    return 0; }
	else if(obc_enth == 2)
	{ disk->bc->obc_enth = FIXED_ENTHALPY_GRADIENT; return 0; }
	else { return -1; }
}

int get_outer_enthalpy_boundary_enthalpy(double *obc_enth_val)
{ *obc_enth_val = disk->bc->obc_enth_val; return 0; }

int set_outer_enthalpy_boundary_enthalpy(double  obc_enth_val)
{ disk->bc->obc_enth_val = obc_enth_val;  return 0; }

int get_outer_enthalpy_boundary_enthalpy_gradient
	(double *obc_enth_val)
{ *obc_enth_val = disk->bc->obc_enth_val; return 0; }

int set_outer_enthalpy_boundary_enthalpy_gradient
	(double  obc_enth_val)
{ disk->bc->obc_enth_val = obc_enth_val;  return 0; }

int get_outer_boundary_function(int *obc_func) 
{ *obc_func = disk->bc->obc_func; return 0; }

int set_outer_boundary_function(int  obc_func) 
{ disk->bc->obc_func = obc_func;  return 0; }


int get_number_of_cells(int *n) { *n = disk->n; return 0; }


int get_number_of_user_parameters(int *n) {
	*n = disk->number_of_parameters;
	return 0;
}

int set_number_of_user_parameters(int n) { 
	if (n > 0) {
		delete disk->params;
		disk->number_of_parameters = n;	
		disk->params = new double[n];
		for (int i = 0; i < n; i++) { disk->params[i] = 0.; }
		return 0; 
	}
	else { return -1; }
}


int get_parameter(int i, double *param) { 
	if (i < 0 || i >= disk->number_of_parameters) {
		return -1;
	}
	else {
		*param = disk->params[i]; 
		return 0; 
	}
}

int set_parameter(int i, double param) { 
	if (i < 0 || i >= disk->number_of_parameters) {
		return -1;
	}
	else {
		disk->params[i] = param;
		return 0; 
	}
}


int get_time(double *time) 
{ *time = disk->current_time; return 0; }


int get_inner_boundary_mass_out(double *mBndOut) {
	*mBndOut = disk->out->mBndOut[0];
	return 0;
}

int get_outer_boundary_mass_out(double *mBndOut) {
	*mBndOut = disk->out->mBndOut[1];
	return 0;
}


int get_inner_boundary_energy_out(double *eBndOut) {
	*eBndOut = disk->out->eBndOut[0];
	return 0;
}

int get_outer_boundary_energy_out(double *eBndOut) {
	*eBndOut = disk->out->eBndOut[1];
	return 0;
}


int get_mass_source_out(int *i, double *mSrcOut, 
						int number_of_points) {

	int i0;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			mSrcOut[j] = disk->out->mSrcOut[i0];
		}
		else { return -1; }
	}

	return 0;
}

int get_energy_source_out(int *i, double *eSrcOut, 
						int number_of_points) {

	int i0;

	for (int j = 0; j < number_of_points; j++) {
		i0 = i[j];

		if (i0 < disk->n) {
			eSrcOut[j] = disk->out->eSrcOut[i0];
		}
		else { return -1; }
	}

	return 0;
}


int get_dtStart(double *dtStart) { 
	*dtStart = disk->ctrl->dtStart;
	return 0;
}

int set_dtStart(double  dtStart) {
	if (dtStart > 0.) { disk->ctrl->dtStart = dtStart; return 0; }
	else { return -1; }
}

int get_dtMin(double *dtMin) {
	*dtMin = disk->ctrl->dtMin;
	return 0;
}

int set_dtMin(double  dtMin) {
	if (dtMin > 0.) { disk->ctrl->dtMin = dtMin; return 0; }
	else { return -1; }
}

int get_dtTol(double *dtTol) {
	*dtTol = disk->ctrl->dtTol;
	return 0;
}

int set_dtTol(double  dtTol) {
	if (dtTol > 0.) { disk->ctrl->dtTol = dtTol; return 0; }
	else { return -1; }
}

int get_errTol(double *errTol) {
	*errTol = disk->ctrl->errTol;
	return 0;
}

int set_errTol(double  errTol) {
	if (errTol > 0.) { disk->ctrl->errTol = errTol; return 0; }
	else { return -1; }
}

int get_maxDtIncrease(double *maxDtIncrease) {
	*maxDtIncrease = disk->ctrl->maxDtIncrease;
	return 0;
}

int set_maxDtIncrease(double  maxDtIncrease) {
	if (maxDtIncrease > 1.) { 
		disk->ctrl->maxDtIncrease = maxDtIncrease;
		return 0;
	}
	else {
		return -1;
	}
}

int get_maxIter(int *maxIter) {
	*maxIter = disk->ctrl->maxIter;
	return 0;
}

int set_maxIter(int  maxIter) {
	if (maxIter > 0) { disk->ctrl->maxIter = maxIter; return 0; }
	else { return -1; }
}

int get_interpOrder(int *interpOrder) {
	*interpOrder = disk->ctrl->interpOrder;
	return 0;
}

int set_interpOrder(int  interpOrder) {
	if (interpOrder > 0 and interpOrder < 4) {
		disk->ctrl->interpOrder = interpOrder;
		return 0;
	}
	else { return -1; }
}

int get_maxStep(int *maxStep) {
	*maxStep = disk->ctrl->maxStep;
	return 0;
}

int set_maxStep(int  maxStep) {
	if (maxStep != 0) { disk->ctrl->maxStep = maxStep; return 0; }
	else { return -1; }
}

int get_useBE(int *useBE) {
	*useBE = disk->ctrl->useBE;
	return 0;
}

int set_useBE(int  useBE) {
	disk->ctrl->useBE = useBE;
	return 0;
}

int get_verbosity(int *verbosity) {
	*verbosity = disk->ctrl->verbosity;
	return 0;
}

int set_verbosity(int  verbosity) {
	if (verbosity > -1 and verbosity < 4) {
		disk->ctrl->verbosity = verbosity;
		return 0;
	}
	else {
		return -1;
	}
}

int get_PreTimestep(int *PreTimestep) {
	*PreTimestep = disk->ctrl->preTimestep_func;
	return 0;
}

int set_PreTimestep(int  PreTimestep) {
	disk->ctrl->preTimestep_func = PreTimestep;
	return 0;
}

int get_PostTimestep(int *PostTimestep) {
	*PostTimestep = disk->ctrl->postTimestep_func;
	return 0;
}

int set_PostTimestep(int  PostTimestep) {
	disk->ctrl->postTimestep_func = PostTimestep;
	return 0;
}

int get_nFail(int *nFail) {
	*nFail = *(disk->out->nFail);
	return 0;
}


int   commit_parameters() { return 0; }

int recommit_parameters() { return 0; }

int cleanup_code() { 
  delete disk; 

  //fprintf(stderr,"Hello from cleanup_code");

  //exit(529);
  //abort();
  //raise(SIGINT);

  return 0; }