#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/***********************************************************************************/
//  Silicon (diamond structure) with Stillinger Weber potential
//
// REFs.:
// Takanobu Watanabe, "Molecular Dynamics and Interatomic Potentials",
// p.102, MORIKITA Syuppann, Tokyo (2023), in Japanese.
// 
// F. H. Stillinger and T. A. Weber, "Computer simulation of local
// order in condensed phases of silicon", Phys. Rev. B., Vol. 31, No. 8, 5262 (1985)
//
/***********************************************************************************/
void get_control_param(void)
{
  ctl.calc_max = 400000;  /* maximum MD time step                     */
  ctl.delta_time_fs = 0.38317;   /* 1 [fs] = 0.001 [ps] = 1.0 X 10^(-15) [s] */

  ctl.temp = 100.0;       /* Temperature setting [K] */
  ctl.t_control_step = 5;
  ctl.set_press_GPa_X = 0.01;  /* Pressure setting  [GPa] */
  ctl.set_press_GPa_Y = 0.01; 
  ctl.set_press_GPa_Z = 0.01; 
  ctl.p_control_step = 20;
  sys.pres_scalling_factor = 0.5;

  sys.Ax = 5.42;   /* initial lattice constant [Si with diamond structure] */
  sys.Ay = 5.42;
  sys.Az = 5.42;
  sys.nx = 5;	   /* number of unit cells in x-direction */
  sys.ny = 5;
  sys.nz = 5;

  ctl.natoms_in_unit_cell = 8;     /* number of atoms in unit cell          */
  ctl.natoms_in_mol_unit  = 2;     /* number of atoms in primitive Mol unit */
  ctl.kinds_of_ions       = 1;     /* Si */

  /* alpha setting for EWALD calculation */
  sys.a1   = 0.0;   /* to SKIP the EWALD calculation */
  sys.hm = 20;      /* has no meaning for SW-Si */

}

/*****
  identify_ion:
  every atom is "Si" in Si crystal with diamond structure 
*****/
void	identify_ion(void)
{
  short	x, y, z, i=0;

  for(z=0; z<sys.nz*2; z++)
    for(y=0; y<sys.ny*2; y++)
      for(x=0; x<sys.nx*2; x++){
	sys.ion[i] = 0;
	i++;
      }
}

/*****
  void set_loc(void)
  set up a diamond structure for Si 
*****/
void	set_loc(void)
{
  short	x, y, z, i=0;
  short nx, ny, nz;
  double nnx=0.0, nny=0.0, nnz=0.0;
  double nx3=0.0, ny3=0.0, nz3=0.0;
  double shiftx, shifty, shiftz;

  nx = sys.nx * 2;
  ny = sys.ny * 2;
  nz = sys.nz * 2;

/**** 
  X display shift 
****/
  shiftx = (sys.Ax)/2.0;
  shifty = (sys.Ay)/2.0;
  shiftz = (sys.Az)/2.0;

/****
  fcc ion
****/

  for(z=0; z<nz; z++)
    for(y=0; y<ny; y++)
      for(x=0; x<nx; x++)
	{
	  if((x+y+z)%2==0){
	    nnx = (double)(x)/2.0;
	    nny = (double)(y)/2.0;
	    nnz = (double)(z)/2.0;
	    sys.rx[i] = sys.Ax * nnx +shiftx;
	    sys.ry[i] = sys.Ay * nny +shifty;
	    sys.rz[i] = sys.Az * nnz +shiftz;
	    i++;

	    nx3 = nnx + 0.25;
	    ny3 = nny + 0.25;
	    nz3 = nnz + 0.25;
	    sys.rx[i] = sys.Ax * nx3 +shiftx;
	    sys.ry[i] = sys.Ay * ny3 +shifty;
	    sys.rz[i] = sys.Az * nz3 +shiftz;
	    i++;
	  }
	}
}

void set_potential(void)
{
  double m_Si;

  m_Si = 28.0855 /(6.0221367 * 1.6605402 / 10.0); /* [m_u] */

  ion.m[0] = m_Si;

  /* ------------ S.W Potential Set --------------*/
  sw.ep = (50.0* 0.04337) * 1.60217733e-1;  /* [eV -> energy] */
  sw.sig = 2.0951;                 /* [A] */
  sw.A = 7.049556277*sw.ep;        /* [energy]       */
  sw.B = 0.6022245584*(sw.sig*sw.sig*sw.sig*sw.sig);
  sw.a = 1.8*sw.sig;               /* cut-off radius */
  sw.lam = 21.0*sw.ep;
  sw.gam = 1.2*sw.sig;
  sw.beta = -2.0*sw.ep*sw.lam;
  /* ---------------------------------------------*/

}

void mk_table(void)
{
  /*  Silicon (diamond structure) with Stillinger Weber potential  */
  /*  this function do nothing.  see "real.c".                     */
}


void   md_xyz(void)
{
  int i;

  fprintf(fpmdxyz,"%6d \n",sys.N);
  fprintf(fpmdxyz,"Lattice=\"%3.6f 0.0 0.0 ",sys.Lx);
  fprintf(fpmdxyz,"0.0 %3.6f 0.0 ",sys.Ly);
  fprintf(fpmdxyz,"0.0 0.0 %3.6f\" ",sys.Lz);
  fprintf(fpmdxyz,"Properties=species:S:1:pos:R:3 %6d\n",sys.step);    
  for(i=0;i<sys.N;i++) { /* [A] unit */
    fprintf(fpmdxyz,"Si ");
  }
    fprintf(fpmdxyz," %3.6f   %3.6f   %3.6f \n", sys.rx[i],sys.ry[i],sys.rz[i]);
}

