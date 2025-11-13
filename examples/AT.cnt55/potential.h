/*********************************************************************/
/*  Carbon (diamond structure) with J.Tersoff potential              */
/*********************************************************************/
FILE *ftube;

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

/*------------------- declaration for the structures ----------------------*/
  potential_at_set at;
