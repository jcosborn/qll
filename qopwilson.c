#include "common.h"
#include <qop.h>
#include <qdp_fn.h>
#include <qdp_dn.h>

void
qopWilsonDslashInit(int *argc, char ***argv)
{
  QDP_initialize(argc, argv);
  QDP_verbose(0);
  QDP_profcontrol(0);
}

void
qopWilsonDslashFini(void)
{
  QDP_finalize();
}

void
qopWilsonDslashSetup(Layout *l)
{
  int nd = l->nDim;
  int *ls = l->physGeom;
  //printf("%i nd: %i\t", myrank, nd);
  //for(int i=0; i<nd; i++) printf(" %i", ls[i]);
  //printf("\n");
  //QDP_finalize();
  //exit(1);
  QDP_set_latsize(nd, ls);
  //QMP_barrier();
  //TRACE_ALL;
  QDP_create_layout();
  //TRACE_ALL;
  //QMP_barrier();
  QOP_layout_t qoplayout = QOP_LAYOUT_ZERO;
  qoplayout.latdim = nd;
  qoplayout.latsize = (int *) malloc(nd*sizeof(int));
  for(int i=0; i<nd; i++) {
    qoplayout.latsize[i] = ls[i];
  }
  qoplayout.machdim = -1;
  QOP_init(&qoplayout);
}

void
toQDP_F(real *xx, QDP_Lattice *lat, real *yy, Layout *l, int nelem, int swap)
{
  int nd = l->nDim;
  int nl = l->nSitesInner;
  int ysize = nelem * nl;
  //double n2=0;
  int nsites = l->nSites;
  if(nsites!=QDP_sites_on_node_L(lat)) {
    printf("%s: nsites(%i) != QDP_sites_on_node_L(lat)(%i)\n",
	   __func__, nsites, QDP_sites_on_node_L(lat));
    QDP_abort(-1);
  }
  for(int i=0; i<nsites; i++) { // QDP sites
    int x[nd];
    LayoutIndex li;
    QDP_get_coords_L(lat, x, QDP_this_node, i);
    layoutIndex(l, &li, x);
    if(li.rank==l->myrank) {
      int oi = li.index / l->nSitesInner;
      int ii = li.index % l->nSitesInner;
      real *yi = yy + (ysize*oi + ii);
      for(int e=0; e<nelem; e++) {
	int e2 = e/2;
	int ei = e2 % swap;
	int eo = e2 / swap;
	int es = ei*(nelem/swap) + eo*2 + (e&1);
	xx[i*nelem+es] = yi[nl*e];
	//n2 += er*er + ei*ei;
      }
    } else {
      printf("unpack: site on wrong node!\n");
      QDP_abort(-1);
    }
  }
  //printf("unpack2: %g\n", n2);
}

void
unpack(Layout *l, int nelem, QDP_N_ColorVector *xx, real *yy)
{
  QDP_Lattice *lat = QDP_N_get_lattice_V(xx);
  real *x = QDP_N_expose_V(xx);
  toQDP_F(x, lat, yy, l, nelem, 1);
  QDP_N_reset_V(xx);
}

void
unpackH(Layout *l, QDP_HalfFermion *xx, real *y)
{
  QDP_Lattice *lat = QDP_get_lattice_H(xx);
  real *x = QDP_expose_H(xx);
  toQDP_F(x, lat, y, l, 12, 3);
  QDP_reset_H(xx);
}

void
unpackD(Layout *l, QDP_DiracFermion *xx, real *y)
{
  QDP_Lattice *lat = QDP_get_lattice_D(xx);
  real *x = QDP_expose_D(xx);
  toQDP_F(x, lat, y, l, 24, 3);
  QDP_reset_D(xx);
}

void
unpackM(Layout *l, QDP_ColorMatrix *m, real *y)
{
  QDP_Lattice *lat = QDP_get_lattice_M(m);
  QDP_ColorMatrix *xx = QDP_create_M_L(lat);
  real *x = QDP_expose_M(xx);
  int nelem = 2*QDP_Nc*QDP_Nc;
  toQDP_F(x, lat, y, l, nelem, 1);
  QDP_reset_M(xx);
  //QDP_M_eq_M(m, xx, QDP_all_L(lat));
  QDP_M_eq_transpose_M(m, xx, QDP_all_L(lat));
  QDP_destroy_M(xx);
}

void
fromQDP_F(real *yy, Layout *l, real *xx, QDP_Lattice *lat, int nelem, int swap)
{
  int nd = l->nDim;
  int nl = l->nSitesInner;
  int ysize = nelem * nl;
  //double n2=0;
  int nsites = l->nSites;
  if(nsites!=QDP_sites_on_node_L(lat)) {
    printf("%s: nsites(%i) != QDP_sites_on_node_L(lat)(%i)\n",
	   __func__, nsites, QDP_sites_on_node_L(lat));
    QDP_abort(-1);
  }
  for(int j=0; j<nsites; j++) { // qll sites
    int x[nd];
    LayoutIndex li;
    li.rank = myrank;
    li.index = j;
    layoutCoord(l, x, &li);
    int r = QDP_node_number_L(lat, x);
    if(r==QDP_this_node) {
      int i = QDP_index_L(lat, x);
      int oi = j / l->nSitesInner;
      int ii = j % l->nSitesInner;
      real *yi = yy + (ysize*oi + ii);
      for(int e=0; e<nelem; e++) {
	int e2 = e/2;
	int ei = e2 % swap;
	int eo = e2 / swap;
	int es = ei*(nelem/swap) + eo*2 + (e&1);
	yi[nl*e] = xx[i*nelem+es];
	//n2 += er*er + ei*ei;
      }
    } else {
      printf("unpack: site on wrong node!\n");
      QDP_abort(-1);
    }
  }
  //printf("unpack2: %g\n", n2);
}

void
pack(Layout *l, int nelem, real *xx, QDP_N_ColorVector *yy)
{
  QDP_Lattice *lat = QDP_N_get_lattice_V(yy);
  real *y = QDP_N_expose_V(yy);
  fromQDP_F(xx, l, y, lat, nelem, 1);
  QDP_N_reset_V(yy);
}

void
packH(Layout *l, real *x, QDP_HalfFermion *yy)
{
  QDP_Lattice *lat = QDP_get_lattice_H(yy);
  real *y = QDP_expose_H(yy);
  fromQDP_F(x, l, y, lat, 12, 3);
  QDP_reset_H(yy);
}

void
packD(Layout *l, real *x, QDP_DiracFermion *yy)
{
  QDP_Lattice *lat = QDP_get_lattice_D(yy);
  real *y = QDP_expose_D(yy);
  fromQDP_F(x, l, y, lat, 24, 3);
  QDP_reset_D(yy);
}

void
packM(Layout *l, real *y, QDP_ColorMatrix *m)
{
  QDP_Lattice *lat = QDP_get_lattice_M(m);
  QDP_ColorMatrix *xx = QDP_create_M_L(lat);
  //QDP_M_eq_M(xx, m, QDP_all_L(lat));
  QDP_M_eq_transpose_M(xx, m, QDP_all_L(lat));
  real *x = QDP_expose_M(xx);
  int nelem = 2*QDP_Nc*QDP_Nc;
  fromQDP_F(y, l, x, lat, nelem, 1);
  QDP_reset_M(xx);
  QDP_destroy_M(xx);
}

void
qopShift(Layout *l, real *xx, real *yy, int nelem, int dir, int len)
{
  int nc = nelem/2;
  QDP_N_ColorVector *x = QDP_N_create_V(nc);
  QDP_N_ColorVector *y = QDP_N_create_V(nc);
  unpack(l, nelem, x, yy);
  int mu = abs(dir) - 1;
  QDP_ShiftDir fb = (dir>0) ? QDP_forward : QDP_backward;
  for(int i=0; i<len; i++) {
    QDP_N_ColorVector *t = x; x = y; y = t;
    QDP_N_V_eq_sV(x, y, QDP_neighbor[mu], fb, QDP_all);
  }
  pack(l, nelem, xx, x);
  QDP_N_destroy_V(x);
  QDP_N_destroy_V(y);
}

void
qopSproj(Layout *l, real *x, real *y, int dir, int sign)
{
  QDP_DiracFermion *in = QDP_create_D();
  QDP_HalfFermion *out = QDP_create_H();
  unpackD(l, in, y);
  QDP_H_eq_spproj_D(out, in, dir, sign, QDP_all);
  packH(l, x, out);
  QDP_destroy_D(in);
  QDP_destroy_H(out);
}

void
qopSrecon(Layout *l, real *x, real *y, int dir, int sign)
{
  QDP_HalfFermion *in = QDP_create_H();
  QDP_DiracFermion *out = QDP_create_D();
  unpackH(l, in, y);
  QDP_D_eq_sprecon_H(out, in, dir, sign, QDP_all);
  packD(l, x, out);
  QDP_destroy_H(in);
  QDP_destroy_D(out);
}

void
qopWilsonDslash(Layout *l, real *x, real *u[8], real mass, int sign,
		real *y, char *sub)
{
  QDP_ColorMatrix *qu[4];
  QDP_DiracFermion *out, *in;
  in = QDP_create_D();
  out = QDP_create_D();
  unpackD(l, in, y);
  unpackD(l, out, x);
  for(int i=0; i<4; i++) {
    qu[i] = QDP_create_M();
    unpackM(l, qu[i], u[2*i]);
    QLA_Real two = 2;
    QDP_M_eq_r_times_M(qu[i], &two, qu[i], QDP_all);
  }
  QOP_FermionLinksWilson *fla;
  fla = QOP_wilson_create_L_from_qdp(qu, NULL);
  QOP_evenodd_t eoOut=QOP_EVENODD, eoIn=QOP_EVENODD;
  if(sub[0]=='e') {
    eoOut = QOP_EVEN;
    eoIn = QOP_ODD;
  }
  if(sub[0]=='o') {
    eoOut = QOP_ODD;
    eoIn = QOP_EVEN;
  }
  real kappa = 0.5/(4+mass);
  QOP_wilson_dslash_qdp(NULL, fla, kappa, sign, out, in, eoOut, eoIn);
  QLA_Real n2;
  QDP_r_eq_norm2_D(&n2, out, QDP_all);
  printf0("out2: %g\n", n2);
  packD(l, x, out);
  QDP_destroy_D(in);
  QDP_destroy_D(out);
  for(int i=0; i<4; i++) {
    QDP_destroy_M(qu[i]);
  }
}

void
qopWilsonSolveMulti(Layout *l, real *x[], real *u[8], double masses[],
		    real *y, int nmasses, double rsq, char *sub)
{
  QDP_ColorMatrix *qu[4];
  QDP_DiracFermion *out[nmasses], *in, **outp;
  outp = out;
  in = QDP_create_D();
  unpackD(l, in, y);
  for(int i=0; i<nmasses; i++) {
    out[i] = QDP_create_D();
    unpackD(l, out[i], x[i]);
    QDP_D_eq_zero(out[i], QDP_even);
  }
  for(int i=0; i<4; i++) {
    qu[i] = QDP_create_M();
    unpackM(l, qu[i], u[2*i]);
    QLA_Real two = 2;
    QDP_M_eq_r_times_M(qu[i], &two, qu[i], QDP_all);
  }
  QOP_FermionLinksWilson *fla;
  fla = QOP_wilson_create_L_from_qdp(qu, NULL);
#if 0
  QOP_evenodd_t eo = QOP_EVENODD;
  if(sub[0]=='e') {
    eo = QOP_EVEN;
  }
  if(sub[0]=='o') {
    eo = QOP_ODD;
  }
#endif
  QOP_evenodd_t eo = QOP_EVEN;
  QOP_info_t info = QOP_INFO_ZERO;
  QOP_invert_arg_t inv_arg = QOP_INVERT_ARG_DEFAULT;
  inv_arg.max_iter = 1000;
  inv_arg.restart = 500;
  inv_arg.max_restarts = 5;
  inv_arg.evenodd = eo;
  inv_arg.mixed_rsq = 0;
  QOP_resid_arg_t res_arg = QOP_RESID_ARG_DEFAULT;
  res_arg.rsqmin = rsq;
  QOP_resid_arg_t *ra[nmasses];
  QOP_resid_arg_t **rap = ra;
  real mf[nmasses], *mfp;
  mfp = mf;
  for(int i=0; i<nmasses; i++) {
    ra[i] = &res_arg;
    mf[i] = masses[i];
  }

  //QOP_verbose(3);
  QOP_wilson_invert_multi_qdp(&info, fla, &inv_arg, &rap,
			      &mfp, &nmasses, &outp, &in, 1);
  //QLA_Real n2;
  //QDP_r_eq_norm2_D(&n2, (QDP_DiracFermion*)out, QDP_all);
  printf0("QOP its: %i\n", res_arg.final_iter);
  QDP_destroy_D(in);
  for(int i=0; i<nmasses; i++) {
    packD(l, x[i], out[i]);
    QDP_destroy_D(out[i]);
  }
  for(int i=0; i<4; i++) {
    QDP_destroy_M(qu[i]);
  }
}

void
qopWilsonSolve(Layout *l, real *x, real *u[8], real mass, real *y,
	       double rsq, char *sub)
{
  QDP_ColorMatrix *qu[4];
  QDP_DiracFermion *out, *in;
  in = QDP_create_D();
  out = QDP_create_D();
  unpackD(l, in, y);
  unpackD(l, out, x);
  for(int i=0; i<4; i++) {
    qu[i] = QDP_create_M();
    unpackM(l, qu[i], u[2*i]);
    QLA_Real two = 2;
    QDP_M_eq_r_times_M(qu[i], &two, qu[i], QDP_all);
  }
  QOP_FermionLinksWilson *fla;
  fla = QOP_wilson_create_L_from_qdp(qu, NULL);
  QOP_evenodd_t eo=QOP_EVENODD;
  if(sub[0]=='e') {
    eo = QOP_EVEN;
  }
  if(sub[0]=='o') {
    eo = QOP_ODD;
  }
  QOP_info_t info = QOP_INFO_ZERO;
  QOP_invert_arg_t inv_arg = QOP_INVERT_ARG_DEFAULT;
  QOP_resid_arg_t res_arg = QOP_RESID_ARG_DEFAULT;
  res_arg.rsqmin = rsq;
  inv_arg.max_iter = 1000;
  inv_arg.restart = 500;
  inv_arg.max_restarts = 5;
  inv_arg.evenodd = eo;
  inv_arg.mixed_rsq = 0;

  QDP_D_eq_zero(out, QDP_even);
  //QOP_verbose(3);
  QOP_wilson_invert_qdp(&info, fla, &inv_arg, &res_arg, mass, out, in);
  //QLA_Real n2;
  //QDP_r_eq_norm2_D(&n2, (QDP_DiracFermion*)out, QDP_all);
  printf0("QOP its: %i\n", res_arg.final_iter);
  packD(l, x, out);
  QDP_destroy_D(in);
  QDP_destroy_D(out);
  for(int i=0; i<4; i++) {
    QDP_destroy_M(qu[i]);
  }
}
