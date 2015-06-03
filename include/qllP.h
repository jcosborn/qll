typedef struct {
  Layout *l;
  real *f;
  int nelem;
  //int nc;
} PR(Field);

void PR(fieldNew)(PR(Field) *f, Layout *l, int nelem);
void PR(fieldFree)(PR(Field) *f);
void PR(fieldSet)(PR(Field) *f, real *r, Layout *l, int nelem/*, int nc*/);
void PR(fieldZero)(PR(Field) *f, Subset *sub);
double PR(fieldSumElem)(PR(Field) *f, Subset *sub);
double PR(fieldNorm2)(PR(Field) *f, Subset *sub);

typedef struct {
  real *u;
  ShiftIndices si;
  int nelem;
  int nc;
  int td;
  ShiftBuf sb;
  Subset sub;
  int backward;
  real *dest;
  real *src;
} PR(Transporter);

void PR(transporterNew)(PR(Transporter) *t, real *u, ShiftIndices *si,
			Subset *sub, int nelem, int nc);
void PR(transporterFree)(PR(Transporter) *t);
void PR(transporterStartT)(PR(Transporter) *t, real *dest, real *src,
			   int td, ThreadCtx *TC);
void PR(transporterLocalT)(PR(Transporter) *t, ThreadCtx *TC);
void PR(transporterWaitT)(PR(Transporter) *t, ThreadCtx *TC);
void PR(transporterDoT)(PR(Transporter) *t, real *dest, real *src, int td,
			ThreadCtx *TC);
void PR(transporterDo)(PR(Transporter) *t, real *dest, real *src, int td);

void PR(gaugeDoublestore)(PR(Field) *ua[], PR(Field) *u[], int dirs[], int lens[], int n);
void PR(gaugePlaq)(double plaqs[], PR(Field) *u[]);
void PR(smearFat7)(PR(Field) *fl[], PR(Field) *ll[], Fat7Coeffs *coef,
		   PR(Field) *u[], PR(Field) *uLong[]);

typedef struct {
  double massFac;
  PR(Transporter) *t;
  int n;
  int fb;
} PR(StaggeredLinks);

typedef struct {
  PR(StaggeredLinks) sle, slo;
  double mass;
  real *v;
  PR(Cgls) cgls;
  int its;
  int noscale;
  double flops;
  double secs;
  Subset sube, subo;
} PR(StaggeredSolver);

void PR(stagDslashSetup)(Layout *layout, int nd, int ls[], int rg[]);
void PR(createStaggeredLinks)(PR(StaggeredLinks) *sl, Layout *l,
			      real *u[], real *u3[], char *sub);
void PR(dslash2)(PR(StaggeredLinks) *d, real *x, real mass, real *y);

void PR(createStaggeredSolver)(PR(StaggeredSolver) *ss, Layout *l,
			       real *u[], real *u3[], char *sub);
void PR(freeStaggeredSolver)(PR(StaggeredSolver) *ss);
int PR(solve)(real *x, real *u[8], real *u3[8], real mass, real *y,
	      double rsq, char *sub);
int PR(solveMulti)(real *x[], real *u[8], real *u3[8], double masses[],
		   int nmasses, real *y, double rsq, char *sub);
int PR(solve2)(PR(StaggeredSolver) *ss, real *x, real mass, real *y,
	       double rsq, int maxits, char *sub);
int PR(solveMulti2)(PR(StaggeredSolver) *ss, real *x[], double masses[],
		    int nmasses, real *y, double rsq, int maxits, char *sub);

typedef struct {
  double massFac;
  PR(Transporter) *t;
  int n;
  int fb;
} PR(WilsonLinks);

typedef struct {
  PR(WilsonLinks) sle, slo;
  double mass;
  real *v;
  PR(Cgls) cgls;
  int its;
  int noscale;
  double flops;
  double secs;
  Subset sube, subo;
} PR(WilsonSolver);

void PR(createWilsonLinks)(PR(WilsonLinks) *sl, Layout *l,
			      real *u[], char *sub);
void PR(wilsonDslash)(PR(WilsonLinks) *d, real *x, real mass, int sign, real *y);

void PR(createWilsonSolver)(PR(WilsonSolver) *ss, Layout *l,
			       real *u[], char *sub);
void PR(freeWilsonSolver)(PR(WilsonSolver) *ss);
int PR(wilsonSolve)(PR(WilsonSolver) *ss, real *x, real mass, real *y,
		    double rsq, int maxits, char *sub);
int PR(wilsonSolveMulti)(PR(WilsonSolver) *ss, real *x[], double masses[],
		    int nmasses, real *y, double rsq, int maxits, char *sub);

