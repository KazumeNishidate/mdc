#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

//  Carbon nanotube with Abell Tersoff potential

/*****
*  purpose: system initialisation
*  all valiables are specified in a MD calculational unit.
*  see also [manual/units.IEMD]
*****/
void get_control_param(void)
{
  char dumy[1000];
  double xx;
  int num=0;

  ctl.calc_max = 500000;        /* maximum MD time step                     */
  ctl.delta_time_fs = 0.38317;  /* 1 [fs] = 0.001 [ps] = 1.0 X 10^(-15) [s] */

  ctl.temp = 100.0;           /* Temperature setting [K]                    */
  ctl.t_control_step = 5;     /* scale at every ctl.t_control_step steps    */
  ctl.set_press_GPa_X = 0.01; /* Pressure setting  [GPa]                    */
  ctl.set_press_GPa_Y = 0.01; 
  ctl.set_press_GPa_Z = 0.01; 
  ctl.p_control_step  =  20;
  sys.pres_scalling_factor = 0.5;

  if((ftube = fopen("./cnt.vasp","r"))==NULL){
    printf("cannot open cnt.vasp. Abort\n");
    exit(EXIT_FAILURE);
  }

  fgets(dumy, 100,ftube); // read the comment
  fgets(dumy, 100,ftube); // read the scaling factor
  fscanf(ftube,"%lf %lf %lf \n",&sys.Ax, &xx, &xx);
  fscanf(ftube,"%lf %lf %lf \n",&xx, &sys.Ay, &xx);
  fscanf(ftube,"%lf %lf %lf \n",&xx, &xx, &sys.Az);
  fgets(dumy, 200,ftube); 
  fscanf(ftube,"%d  \n",&num); // read the number of atoms

  fgets(dumy, 200,ftube);
  ctl.natoms_in_unit_cell = num;
  ctl.natoms_in_mol_unit  = num;
  sys.N = num;

  sys.nx = 1;	   /* number of unit cells in x-direction                 */
  sys.ny = 1;      /* note: currently only nx=ny=nz=1 of C60 is supported */
  sys.nz = 1;

  ctl.kinds_of_ions       = 1;     /* C */

  /* alpha setting for EWALD calculation */
  sys.a1 = 0.0;     /* to SKIP the EWALD calculation */
  sys.hm = 20;      /* has no meaning for AT-C */
  sys.radius = 3.0; /* [A] AT potential set */  
}

/*****
  identify_ion:
  every atom is "C" in C crystal with diamond structure 
*****/
void	identify_ion(void)

{
  short	i=0;

  for(i=0; i<sys.N; i++){
    sys.ion[i] = 0;
  }
}


/*****
  void set_loc(void)
*****/
void	set_loc(void){
  int i;
  double rxx, ryy, rzz;
  
  for(i=0;i<sys.N;i++){
    fscanf(ftube," %lf %lf %lf \n",&rxx, &ryy, &rzz);
    sys.rx[i] = rxx;
    sys.ry[i] = ryy;
    sys.rz[i] = rzz;    
  }
}

void set_potential(void)
{
  double m_C;

  m_C = 12.011 /(6.0221367 * 1.6605402 / 10.0); /* [m_u] */

  ion.m[0] = m_C; // all Carbon

  // vdW set for carbon systems
  //
  // - C. Li and TW. Chou,
  //   Composites Sci. and Techol. 63, 1517-1524 (2003).
  // - H. Hemmatian, M.R. Zamani, JE Jam,
  //   J. Theoret. Appl. Mech. 57, 207 (2019).
  // - L. Battezzati, C. Pisani, and F. Ricca,
  //   J. Chem. Soc. 71, 1629-1639 (1975).

  vdw.sig = 3.4;                     // [angstrom]
  vdw.eps = 0.0556 * 6.94769546e-3;  // [kcal/mol -> en]


  vdw.cut_out = 2.5*vdw.sig; // definition of the inner and the outer cutoff radii
  vdw.cut_in  = 2.8;         //     F_vdW active region {r}:  cut_in < r < cut_out
  
  /* ------------ A.T Potential Set --------------*/
  at.A = 1.3936e3 * 1.60217733e-1;  /* [eV -> en] */
  at.B = 3.4670e2 * 1.60217733e-1;  /* [eV -> en] */
  at.lam = 3.4879;
  at.mu = 2.2119;
  at.beta = 1.5724e-7;
  at.n = 7.2751e-1;
  at.c = 3.8049e4;
  at.d = 4.384;
  at.h = -0.57058;
  at.R = 1.8;
  at.S = 2.1;               /* cut-off radius of AT */
  at.x = 1.0;
  /* ---------------------------------------------*/

}

void mk_table(void)
{
  /*  Carbon (diamond structure) with Abell Tersoff potential      */
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
    fprintf(fpmdxyz,"C %3.6f   %3.6f   %3.6f \n", sys.rx[i],sys.ry[i],sys.rz[i]);
  }
}





