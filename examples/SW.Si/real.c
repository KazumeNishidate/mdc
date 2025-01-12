#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/*****
  void	real_space(void)

  The function for potential and force calculation of Si diamond structure
  with Stillinger Weber potentail set.

*****/
void	real_space(void)
{
  int i, j, k; 
  double rdx, rdy, rdz;
  double x_ij, y_ij, z_ij;
  double x_ik, y_ik, z_ik;
  double x_jk, y_jk, z_jk;
  double r_ij, r_ik, r_jk;

  double r_ij4;

  double for_2, f2_1, f2_2;
  double cs_jik;
  double ex_jik;

  double fa_jik, fb_jik, fd_jik;

  rdx = sys.Lx*0.5;
  rdy = sys.Ly*0.5;
  rdz = sys.Lz*0.5;

  /* initialization */
  sys.pot = 0.0;
  sys.virX = 0.0;  sys.virY = 0.0;  sys.virZ = 0.0;
  
  for(i=0; i<sys.N; i++) {
    sys.fx[i] = 0.0; sys.fy[i] = 0.0; sys.fz[i] = 0.0;
  }
  
  for(i=0; i<sys.N; i++) {
    for(j=i+1; j<sys.N; j++) {
      x_ij = sys.rx[i] - sys.rx[j];
      y_ij = sys.ry[i] - sys.ry[j];
      z_ij = sys.rz[i] - sys.rz[j];
      
      /* cyclic boundary condition */
      if(x_ij > rdx) x_ij -= sys.Lx;
      if(y_ij > rdy) y_ij -= sys.Ly;
      if(z_ij > rdz) z_ij -= sys.Lz;
      if(x_ij < -rdx) x_ij += sys.Lx;
      if(y_ij < -rdy) y_ij += sys.Ly;
      if(z_ij < -rdz) z_ij += sys.Lz;

      /*  ion_i - ion_j distance  */
      r_ij = sqrt(x_ij*x_ij + y_ij*y_ij + z_ij*z_ij);
      
      /* sys.radius = cut-off radius */
      if(r_ij >= sw.a) continue;

      r_ij4 = r_ij * r_ij * r_ij * r_ij;
      
      /* ----------------  SET 2-BODY POTENTIAL ---------------------*/
      sys.pot += sw.ep*sw.A*(sw.B/r_ij4-1.0)*exp(1.0/(r_ij-sw.a));
      
      /* ----------------  SET 2-BODY FORCE ------------------------*/
      f2_1 = -sw.ep*sw.A*exp(1.0/(r_ij-sw.a))/r_ij;
      f2_2 = -(4.0*sw.B/(r_ij4*r_ij)) + (1.0-sw.B/r_ij4)/((r_ij-sw.a)*(r_ij-sw.a));
      for_2 = f2_1 * f2_2;

      sys.fx[i] += (for_2 * x_ij);
      sys.fy[i] += (for_2 * y_ij);
      sys.fz[i] += (for_2 * z_ij);

      sys.fx[j] -= (for_2 * x_ij);
      sys.fy[j] -= (for_2 * y_ij);
      sys.fz[j] -= (for_2 * z_ij);

      /* VIRIAL */
      sys.virX += (for_2 * x_ij) * x_ij; 
      sys.virY += (for_2 * y_ij) * y_ij; 
      sys.virZ += (for_2 * z_ij) * z_ij; 
    }
  }

  /* ------------------  SET 3-BODY ------------------------------*/
  for(i=0;i<sys.N;i++)  {
    for(j=0;j<sys.N;j++)  {
      if(i==j) continue;
      x_ij = sys.rx[i] - sys.rx[j];
      y_ij = sys.ry[i] - sys.ry[j];
      z_ij = sys.rz[i] - sys.rz[j];
      
      /* cyclic boundary condition */
      if(x_ij > rdx) x_ij -= sys.Lx;
      if(y_ij > rdy) y_ij -= sys.Ly;
      if(z_ij > rdz) z_ij -= sys.Lz;
      if(x_ij < -rdx) x_ij += sys.Lx;
      if(y_ij < -rdy) y_ij += sys.Ly;
      if(z_ij < -rdz) z_ij += sys.Lz;

      /*  ion_i - ion_j distance  */
      r_ij = sqrt(x_ij*x_ij + y_ij*y_ij + z_ij*z_ij);
      
      /* sys.radius = cut-off radius */
      if(r_ij >= sw.a) continue;

      for(k=j+1;k<sys.N;k++){
	if(k==i) continue;
	x_ik = sys.rx[i] - sys.rx[k];  
	y_ik = sys.ry[i] - sys.ry[k];  
	z_ik = sys.rz[i] - sys.rz[k];  
	
	/* cyclic boundary condition */
	if(x_ik > rdx) x_ik -= sys.Lx;
	if(y_ik > rdy) y_ik -= sys.Ly;
	if(z_ik > rdz) z_ik -= sys.Lz;
	if(x_ik < -rdx) x_ik += sys.Lx;
	if(y_ik < -rdy) y_ik += sys.Ly;
	if(z_ik < -rdz) z_ik += sys.Lz;
	
	/* ion_i - ion_k distance */
	r_ik = sqrt(x_ik*x_ik + y_ik*y_ik + z_ik*z_ik);
 	
	if(r_ik > sw.a) continue;
	
	x_jk = x_ik - x_ij;
	y_jk = y_ik - y_ij;
	z_jk = z_ik - z_ij;	

	r_jk = sqrt(x_jk*x_jk + y_jk*y_jk + z_jk*z_jk);	
	cs_jik = (x_ij*x_ik + y_ij*y_ik + z_ij*z_ik)/(r_ij * r_ik);
	ex_jik = exp(sw.gam/(r_ij-sw.a)+sw.gam/(r_ik-sw.a));

	sys.pot += sw.lam*ex_jik*(cs_jik+1.0/3.0)*(cs_jik+1.0/3.0);

	fa_jik = -sw.lam*ex_jik
	  *(cs_jik+1.0/3.0)*((cs_jik+1.0/3.0)*sw.gam/(r_ij-sw.a)/(r_ij-sw.a)
			    -2.0*(r_ij-r_ik*cs_jik)/r_ij/r_ik);

	fb_jik = -sw.lam*ex_jik
	  *(cs_jik+1.0/3.0)*((cs_jik+1.0/3.0)*sw.gam/(r_ik-sw.a)/(r_ik-sw.a)
			    -2.0*(r_ik-r_ij*cs_jik)/r_ik/r_ij);

	fd_jik = -sw.lam*exp(sw.gam/(r_ij-sw.a)+sw.gam/(r_ik-sw.a))
	  *(cs_jik+1.0/3.0)*2.0*r_jk/r_ij/r_ik;

	sys.fx[i] += -fa_jik*x_ij/r_ij - fb_jik*x_ik/r_ik;
	sys.fy[i] += -fa_jik*y_ij/r_ij - fb_jik*y_ik/r_ik;
	sys.fz[i] += -fa_jik*z_ij/r_ij - fb_jik*z_ik/r_ik;

	sys.fx[j] +=  fa_jik*x_ij/r_ij - fd_jik*x_jk/r_jk;
	sys.fy[j] +=  fa_jik*y_ij/r_ij - fd_jik*y_jk/r_jk;
	sys.fz[j] +=  fa_jik*z_ij/r_ij - fd_jik*z_jk/r_jk;

	sys.fx[k] +=  fb_jik*x_ik/r_ik + fd_jik*x_jk/r_jk;
	sys.fy[k] +=  fb_jik*y_ik/r_ik + fd_jik*y_jk/r_jk;
	sys.fz[k] +=  fb_jik*z_ik/r_ik + fd_jik*z_jk/r_jk;
	
      }
    }
  }
}

