/*********************************************************************/
/*  Carbon (diamond structure) with J.Tersoff potential              */
/*********************************************************************/

typedef struct{
  double A;
  double B;
  double lam;
  double mu;
  double beta;
  double n;
  double c;
  double d;
  double h;
  double R;
  double S;
  double x;

} potential_at_set;

typedef struct{
  double sig;
  double eps;
  double cut_out;
  double cut_in;
} vdw_set;

// real.c 
void   check_vdw(void);
void   calc_vdw(void);

/*------------------- declaration for the structures ----------------------*/
potential_at_set at;
vdw_set vdw;
