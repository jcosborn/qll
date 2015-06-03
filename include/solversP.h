typedef struct {
  real *r1;
  real *r2;
  real *z;
  real *p;
  real *Ap;
  real *Adr2;
  int nelem;
  int verbose;
  int indent;
} PR(Cgls);

typedef void (PR(Linop))(real *Ax, real *x, LinopInfo *li,
			 void *args, ThreadCtx *TC);

void PR(cglsInit)(PR(Cgls) *cgls, int nelem);
		 void PR(cglsFree)(PR(Cgls) *cgls);
int PR(cgSolve)(PR(Cgls) *cgls, real *x, real *b, PR(Linop) *Aop, void *Aargs,
		PR(Linop) *Mop, void *Margs, double rsq, int itnlim);
int PR(cgmsSolve)(PR(Cgls) *cgls, real *x[], real *b, double shifts[],
		  int nshifts, PR(Linop) *Aop, void *Aargs, PR(Linop) *Mop,
		  void *Margs, double rsqs[], int itnlim);
