#include <common.h>

//M_peq_M(staple, tmat1);
//M_eq_r_times_M(fl[dir], coef1, gf[dir]);
//M_peq_r_times_M(fl[mu], coef, staple);
//M_eq_M_times_M(tmat2, tmat1, ts1);
//M_eq_M_times_Ma(staple, tmat1, ts1);
//M_eq_Ma_times_M(tmat1, gauge[nu], link);

#define M_peq_M(x,y) X_vpeq_r_times_X((x)->f,1,(y)->f,(x)->nelem,ti0,ti1)
#define M_eq_r_times_M(x,r,y) X_veq_r_times_X((x)->f,r,(y)->f,(x)->nelem,ti0,ti1)
#define M_peq_r_times_M(x,r,y) X_vpeq_r_times_X((x)->f,r,(y)->f,(x)->nelem,ti0,ti1)
#define M_eq_M_times_M(x,y,z) do { F_eq_r(x,0); \
    P(X_peq_M_times_X)((x)->f,(y)->f,(z)->f,(x)->nelem,nc,2,ti0,ti1); } while(0)
#define M_eq_M_times_Ma(x,y,z) do { F_eq_r(x,0); \
    P(X_peq_M_times_Xa)((x)->f,(y)->f,(z)->f,(x)->nelem,nc,2,ti0,ti1);}while(0)
#define M_eq_Ma_times_M(x,y,z) do { F_eq_r(x,0); \
    P(X_peq_Ma_times_X)((x)->f,(y)->f,(z)->f,(x)->nelem,nc,2,ti0,ti1); } while(0)

#if 0
static void
dump(real *x, int n, int vl)
{
  for(int i=0; i<n; i++) {
    printf("%12g", x[i*vl]);
    for(int j=1; j<vl; j++) {
      printf(" %12g", x[i*vl+j]);
    }
    printf("\n");
  }
}
#endif

void
gaugeDoublestore(Field *ua[], Field *u[], int dirs[], int lens[], int n)
{
  Transporter t[n];
  char *csub = "all";
  Subset sub;
  layoutSubset(&sub, u[0]->l, csub);
  for(int i=0; i<n; i++) {
    ShiftIndices *si = P(getShift)(dirs[i], lens[i], csub);
    transporterNew(&t[i], NULL, si, &sub, u[i]->nelem, 0);
  }
  beginParallelSubset(sub) {
    for(int i=0; i<n; i++) {
      transporterStartT(&t[i], ua[i]->f, u[i]->f, 0, TC);
    }
    for(int i=0; i<n; i++) {
      transporterLocalT(&t[i], TC);
    }
    for(int i=0; i<n; i++) {
      SUBSET_DEFS2;
      transporterWaitT(&t[i], TC);
      setOffSubC2(sub);
      int nc = floor(0.5+sqrt(0.5*u[i]->nelem));
      Xadj(ua[i]->f, nc, ti0, ti1);
    }
  } endParallelSubset;
  for(int i=0; i<n; i++) {
    transporterFree(&t[i]);
  }
}

static void
setFwdTransporters(Transporter t[], Field *u[], int nelem, int nc, char *csub)
{
  int nd = u[0]->l->nDim;
  Subset sub;
  layoutSubset(&sub, u[0]->l, csub);
  for(int i=0; i<nd; i++) {
    ShiftIndices *si = P(getShift)(i, 1, csub);
    transporterNew(&t[i], u[i]->f, si, &sub, nelem, nc);
  }
}

static void
setBckTransporters(Transporter t[], Field *u[], int nelem, int nc, char *csub)
{
  int nd = u[0]->l->nDim;
  Subset sub;
  layoutSubset(&sub, u[0]->l, csub);
  for(int i=0; i<nd; i++) {
    ShiftIndices *si = P(getShift)(i, -1, csub);
    transporterNew(&t[i], u[i]->f, si, &sub, nelem, nc);
  }
}

#if 0
static void
setTransporters(Transporter t[], Field *u[], int nelem, int nc, char *csub)
{
  int nd = u[0]->l->nDim;
  Subset sub;
  layoutSubset(&sub, u[0]->l, csub);
  Field ua[nd], *uap[nd];
  int dirs[nd], lens[nd];
  for(int i=0; i<nd; i++) {
    ShiftIndices *si = P(getShift)(i, 1, csub);
    transporterNew(&t[2*i], u[i]->f, si, &sub, nelem, nc);
    fieldNew(&ua[i], u[i]->l, u[i]->nelem);
    uap[i] = &ua[i];
    dirs[i] = i;
    lens[i] = -1;
  }
  gaugeDoublestore(uap, u, dirs, lens, nd);
  for(int i=0; i<nd; i++) {
    ShiftIndices *si = P(getShift)(i, -1, csub);
    transporterNew(&t[2*i+1], ua[i].f, si, &sub, nelem, nc);
  }
}
#endif

void
gaugePlaq(double plaqs[], Field *u[])
{
  int nd = u[0]->l->nDim;
  int nc = floor(0.5+sqrt(0.5*u[0]->nelem));
  int nelemM = 2*nc*nc;
  Field f[nd-1][nd];
  Transporter t[nd-1][nd];
  for(int i=0; i<nd-1; i++) {
    for(int j=0; j<nd; j++) {
      fieldNew(&f[i][j], u[0]->l, u[0]->nelem);
    }
    setFwdTransporters(t[i], u, nelemM, nc, "all");
  }
  beginParallelSubset(t[0][0].sub) {
    for(int i=0; i<nd; i++) {
      int k = 0;
      for(int j=0; j<nd; j++) {
	if(j==i) continue;
	transporterStartT(&t[k][i], f[k][i].f, u[j]->f, 2, TC);
	k++;
      }
    }
    for(int i=0; i<nd; i++) {
      int k = 0;
      for(int j=0; j<nd; j++) {
	if(j==i) continue;
	F_eq_r(&f[k][i], 0);
	transporterLocalT(&t[k][i], TC);
	k++;
      }
    }
    SUBSET_DEFS2;
    setOffSubC2(t[0][0].sub);
    int np = (nd*(nd-1))/2;
    double *p, pp[np];
    if(TC->tid==0) p = plaqs;
    else p = pp;
    int l = 0;
    for(int i=0; i<nd; i++) {
      int k = 0;
      for(int j=0; j<nd; j++) {
	if(j==i) continue;
	transporterWaitT(&t[k][i], TC);
	k++;
	if(j<i) { // have reversed elbow
	  //if(u[0]->l->myrank==0) {
	    //dump(f[j][i].f, nelemM, u[0]->l->nSitesInner);
	    //dump(f[j][i].f, u[0]->l->nSitesOuter*nelemM, u[0]->l->nSitesInner);
	    //printf("\n");
	    //dump(f[i-1][j].f, u[0]->l->nSitesOuter*nelemM, u[0]->l->nSitesInner);
	  //}
	  p[l] = re_X_dot_X(f[j][i].f, f[i-1][j].f, nelemM, ti0, ti1);
	  //p[l] = lnorm2_F(&f[j][i]);
	  //p[l] = lnorm2_F(&f[i-1][j]);
	  p[l] /= nc*u[0]->l->physVol;
	  //printf0("%i: %i %g\n", TC->tid, l, p[l]);
	  l++;
	}
      }
    }
    sumDoubleArray(p, np);
  } endParallelSubset;
  for(int i=0; i<nd-1; i++) {
    for(int j=0; j<nd; j++) {
      fieldFree(&f[i][j]);
      transporterFree(&t[i][j]);
    }
  }
}

  //#define printnorm(x) {double n2 = lnorm2_F(x); sumDoubleArray(&n2,1); if(TC->tid==0) { printf0("%s2: %g\n", #x, n2); }}
#define printnorm(x) {double n2 = sumElem_X((x)->f, (x)->nelem, ti0, ti1); sumDoubleArray(&n2,1); if(TC->tid==0) { printf0("%s2: %g\n", #x, n2); }}

static void
genStaple(Field *staple, int mu, int nu,
	  Field *link, double coef,
	  Field **gauge, Field **fl,
	  Field *ts0, Field *ts1,
	  Field *tmat1, Field *tmat2,
	  Transporter tf[], Transporter tb[], ThreadCtx *TC)
{
  int nc = tf[0].nc;
  SUBSET_DEFS2;
  setOffSubC2(tf[0].sub);
  /* Upper staple */
  //if(link!=gauge[mu])
  //QDP_M_eq_sM(ts0, link, QDP_neighbor[nu], QDP_forward, QDP_all);
  //QDP_M_eq_sM(ts1, gauge[nu], QDP_neighbor[mu], QDP_forward, QDP_all);

  if(staple!=NULL) {  /* Save the staple */
    //QDP_M_eq_M_times_Ma(tmat1, ts0, ts1, QDP_all);
    //QDP_M_eq_M_times_M(staple, gauge[nu], tmat1, QDP_all);
    transporterStartT(&tf[nu], tmat1->f, link->f, 2, TC);
    F_eq_r(tmat1, 0);
    transporterLocalT(&tf[nu], TC);
    transporterWaitT(&tf[nu], TC);
    M_eq_M_times_Ma(staple, tmat1, ts1);
  } else {  /* No need to save the staple. Add it to the fatlinks */
    //QDP_M_eq_M_times_Ma(tmat1, ts0, ts1, QDP_all);
    //QDP_M_eq_M_times_M(tmat2, gauge[nu], tmat1, QDP_all);
    //QDP_M_peq_r_times_M(fl[mu], &coef, tmat2, QDP_all);
#if 1
    transporterStartT(&tf[nu], tmat1->f, link->f, 2, TC);
    F_eq_r(tmat1, 0);
    transporterLocalT(&tf[nu], TC);
    transporterWaitT(&tf[nu], TC);
    //printnorm(ts1);
    //printnorm(tmat1);
    M_eq_M_times_Ma(tmat2, tmat1, ts1);
    //printnorm(tmat2);
#else
    transporterStartT(&tf[nu], tmat2->f, link->f, 0, TC);
    transporterLocalT(&tf[nu], TC);
    transporterWaitT(&tf[nu], TC);
    //printnorm(ts1);
    //printnorm(tmat2);
    M_eq_M_times_Ma(tmat1, tmat2, ts1);
    //printnorm(tmat1);
    //printnorm(gauge[nu]);
    M_eq_M_times_M(tmat2, gauge[nu], tmat1);
#endif
    //printnorm(tmat2);
    M_peq_r_times_M(fl[mu], coef, tmat2);
    //printnorm(fl[mu]);
  }

  /* lower staple */
  //QDP_M_eq_sM(ts1, gauge[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
  //QDP_M_eq_Ma_times_M(tmat1, gauge[nu], link, QDP_all);
  //QDP_M_eq_M_times_M(tmat2, tmat1, ts1, QDP_all);
  //QDP_M_eq_sM(ts2, tmat2, QDP_neighbor[nu], QDP_backward, QDP_all);
  M_eq_Ma_times_M(tmat1, gauge[nu], link);
  M_eq_M_times_M(tmat2, tmat1, ts1);
  TBARRIER;
  transporterStartT(&tb[nu], tmat1->f, tmat2->f, 0, TC);
  transporterLocalT(&tb[nu], TC);
  transporterWaitT(&tb[nu], TC);

  if(staple!=NULL) { /* Save the staple */
    //QDP_M_peq_M(staple, ts2, QDP_all);
    //QDP_M_peq_r_times_M(fl[mu], &coef, staple, QDP_all);
    M_peq_M(staple, tmat1);
    M_peq_r_times_M(fl[mu], coef, staple);
  } else { /* No need to save the staple. Add it to the fatlinks */
    //QDP_M_peq_r_times_M(fl[mu], &coef, ts2, QDP_all);
    M_peq_r_times_M(fl[mu], coef, tmat1);
  }
}

void
smearFat7(Field *fl[], Field *ll[], Fat7Coeffs *coef,
	  Field *u[], Field *uLong[])
{
  int nd = u[0]->l->nDim;
  int nc = floor(0.5+sqrt(0.5*u[0]->nelem));
  int nelemM = 2*nc*nc;
  Field staple, tempmat1, t1, t2, tsg[nd][nd];
  Transporter tf[nd], tb[nd];
  //double nflop = 0;
  //double dtime;

  double coef1 = coef->one_link;
  double coef3 = coef->three_staple;
  double coef5 = coef->five_staple;
  double coef7 = coef->seven_staple;
  double coefL = coef->lepage;
  /* to fix up the Lepage term, included by a trick below */
  coef1 -= 6*coefL;
  int have5 = coef5 || coef7 || coefL;
  int have3 = coef3 || have5;

  if(have3 || ll) {
    //nflop = 61632;
    fieldNew(&staple, u[0]->l, u[0]->nelem);
    fieldNew(&tempmat1, u[0]->l, u[0]->nelem);
    if(have3) {
      fieldNew(&t1, u[0]->l, u[0]->nelem);
      fieldNew(&t2, u[0]->l, u[0]->nelem);
      for(int dir=0; dir<nd; dir++) {
        for(int nu=0; nu<nd; nu++) {
          //tsg[dir][nu] = NULL;
          if(dir!=nu) {
            //tsg[dir][nu] = QDP_create_M();
	    fieldNew(&tsg[dir][nu], u[0]->l, u[0]->nelem);
          }
        }
      }
    }
  }
  setFwdTransporters(tf, u, nelemM, nc, "all");
  setBckTransporters(tb, u, nelemM, nc, "all");

  //dtime = -QOP_time();

  beginParallelSubset(tf[0].sub) {
    SUBSET_DEFS2;
    setOffSubC2(tf[0].sub);
    if(have3 || ll) {
      if(have3) {
	for(int dir=0; dir<nd; dir++) {
	  for(int nu=0; nu<nd; nu++) {
	    if(dir!=nu) {
	      //QDP_M_eq_sM(tsg[dir][nu], gf[dir], QDP_neighbor[nu], QDP_forward, QDP_all);
	      transporterStartT(&tf[nu], tsg[dir][nu].f, u[dir]->f, 0, TC);
	      transporterLocalT(&tf[nu], TC);
	      transporterWaitT(&tf[nu], TC);
	    }
	  }
	  TBARRIER;
	}
      }
    }
    for(int dir=0; dir<nd; dir++) {
      //QDP_M_eq_r_times_M(fl[dir], &coef1, gf[dir], QDP_all);
      //printnorm(u[dir]);
      M_eq_r_times_M(fl[dir], coef1, u[dir]);
      //printnorm(fl[dir]);
      if(have3) {
	for(int nu=0; nu<nd; nu++) {
	  if(nu==dir) continue;
	  genStaple(&staple, dir, nu, u[dir], coef3, u, fl,
		    &tsg[dir][nu], &tsg[nu][dir], &t1, &t2, tf, tb, TC);
	  //printnorm(&staple);
          //printnorm(fl[dir]);
	  if(coefL) {
	    TBARRIER;
	    genStaple(NULL, dir, nu, &staple, coefL, u, fl,
		      NULL, &tsg[nu][dir], &t1, &t2, tf, tb, TC);
	    //printnorm(fl[dir]);
	  }
	  if(coef5 || coef7) {
	    for(int rho=0; rho<nd; rho++) {
	      if((rho==dir)||(rho==nu)) continue;
	      TBARRIER;
	      genStaple(&tempmat1, dir, rho, &staple, coef5, u, fl,
			NULL, &tsg[rho][dir], &t1, &t2, tf, tb, TC);
	      //printnorm(&tempmat1);
	      //printnorm(fl[dir]);
	      if(coef7) {
		for(int sig=0; sig<nd; sig++) {
		  if((sig!=dir)&&(sig!=nu)&&(sig!=rho)) {
		    TBARRIER;
		    genStaple(NULL, dir, sig, &tempmat1, coef7,u,fl,
			      NULL, &tsg[sig][dir], &t1, &t2, tf, tb, TC);
		  }
		} /* sig */
	      }
	    } /* rho */
          }
        } /* nu */
      }
    } /* dir */

    /* long links */
    if(ll) {
      double naik = coef->naik;
      for(int dir=0; dir<nd; dir++) {
	//QDP_M_eq_sM(staple,gfLong[dir],QDP_neighbor[dir],QDP_forward,QDP_all);
	//QDP_M_eq_M_times_M(tempmat1, gfLong[dir], staple, QDP_all);
	//QDP_M_eq_sM(staple, tempmat1, QDP_neighbor[dir], QDP_forward, QDP_all);
	//QDP_M_eq_M_times_M(ll[dir], gfLong[dir], staple, QDP_all);
	//QDP_M_eq_r_times_M(ll[dir], &naik, ll[dir], QDP_all);
	tf[dir].u = uLong[dir]->f;
	//xport(tf[dir], tempmat1, uLong[dir]);
	TBARRIER;
	transporterStartT(&tf[dir], tempmat1.f, uLong[dir]->f, 2, TC);
	F_eq_r(&tempmat1, 0);
	transporterLocalT(&tf[dir], TC);
	transporterWaitT(&tf[dir], TC);
	//xport(tf[dir], ll[dir], tempmat1);
	TBARRIER;
	transporterStartT(&tf[dir], ll[dir]->f, tempmat1.f, 2, TC);
	F_eq_r(ll[dir], 0);
	transporterLocalT(&tf[dir], TC);
	transporterWaitT(&tf[dir], TC);
	M_eq_r_times_M(ll[dir], naik, ll[dir]);
      }
    }
  } endParallelSubset;

  if(have3 || ll) {
    fieldFree(&staple);
    fieldFree(&tempmat1);
    if(have3) {
      fieldFree(&t1);
      fieldFree(&t2);
      for(int dir=0; dir<nd; dir++) {
        for(int nu=0; nu<nd; nu++) {
          if(dir!=nu) fieldFree(&tsg[dir][nu]);
        }
      }
    }
  }

  //dtime += QOP_time();
  //node0_printf("LLTIME(Fat): time = %e (Asqtad opt) mflops = %e\n",dtime,
  //       (Real)nflop*volume/(1e6*dtime*numnodes()) );

  //info->final_sec = dtime;
  //info->final_flop = nflop*QDP_sites_on_node;
  //info->status = QOP_SUCCESS;
}
