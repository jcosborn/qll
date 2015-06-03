#include "common.h"
#include <qop.h>
#include <qdp_fn.h>
#include <qdp_dn.h>

void
qopStagDslashInit(int *argc, char ***argv)
{
  QDP_initialize(argc, argv);
  QDP_verbose(0);
  QDP_profcontrol(0);
}

void
qopStagDslashFini(void)
{
  QDP_finalize();
}

void
qopStagDslashSetup(Layout *l)
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
toQDP_F(real *xx, QDP_Lattice *lat, real *yy, Layout *l, int nelem)
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
	xx[i*nelem+e] = yi[nl*e];
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
  toQDP_F(x, lat, yy, l, nelem);
  QDP_N_reset_V(xx);
}

void
unpackV(Layout *l, QDP_ColorVector *xx, real *y)
{
  QDP_Lattice *lat = QDP_get_lattice_V(xx);
  real *x = QDP_expose_V(xx);
  toQDP_F(x, lat, y, l, 6);
  QDP_reset_V(xx);
}

void
unpackM(Layout *l, QDP_ColorMatrix *m, real *y)
{
  QDP_Lattice *lat = QDP_get_lattice_M(m);
  QDP_ColorMatrix *xx = QDP_create_M_L(lat);
  real *x = QDP_expose_M(xx);
  int nelem = 2*QDP_Nc*QDP_Nc;
  toQDP_F(x, lat, y, l, nelem);
  QDP_reset_M(xx);
  //QDP_M_eq_M(m, xx, QDP_all_L(lat));
  QDP_M_eq_transpose_M(m, xx, QDP_all_L(lat));
  QDP_destroy_M(xx);
}

void
fromQDP_F(real *yy, Layout *l, real *xx, QDP_Lattice *lat, int nelem)
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
	yi[nl*e] = xx[i*nelem+e];
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
  fromQDP_F(xx, l, y, lat, nelem);
  QDP_N_reset_V(yy);
}

void
packV(Layout *l, real *x, QDP_ColorVector *yy)
{
  QDP_Lattice *lat = QDP_get_lattice_V(yy);
  real *y = QDP_expose_V(yy);
  fromQDP_F(x, l, y, lat, 6);
  QDP_reset_V(yy);
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
  fromQDP_F(y, l, x, lat, nelem);
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
qopTransportM(Layout *l, real *xx, real *mm, real *yy, int dir, int len)
{
  QDP_ColorMatrix *x = QDP_create_M();
  QDP_ColorMatrix *m = QDP_create_M();
  QDP_ColorMatrix *y = QDP_create_M();
  unpackM(l, x, xx);
  unpackM(l, m, mm);
  unpackM(l, y, yy);
  if(len>0) {
    QDP_M_eq_M_times_sM(x, m, y, QDP_neighbor[dir], QDP_forward, QDP_all);
  } else {
    QDP_M_eq_Ma_times_sM(x, m, y, QDP_neighbor[dir], QDP_backward, QDP_all);
  }
  packM(l, xx, x);
  QDP_destroy_M(x);
  QDP_destroy_M(m);
  QDP_destroy_M(y);
}

void
qopSmear(Field *fl[], Field *ll[], Fat7Coeffs *coef, Field *u[], Field *ul[])
{
  Layout *l = fl[0]->l;
  QDP_ColorMatrix *qu[4], *qul[4];
  for(int i=0; i<4; i++) {
    qu[i] = QDP_create_M();
    unpackM(l, qu[i], u[i]->f);
    if(ul && ul!=u) {
      qul[i] = QDP_create_M();
      unpackM(l, qul[i], ul[i]->f);
    }
  }
  QOP_GaugeField *qg = QOP_create_G_from_qdp(qu);
  QOP_GaugeField *qg2 = qg;
  if(ul && ul!=u) {
    qg2 = QOP_create_G_from_qdp(qul);
  }
  QOP_info_t info;
  QOP_FermionLinksAsqtad *fla;
  fla = QOP_asqtad_create_L_from_G2(&info,(QOP_asqtad_coeffs_t*)coef,qg,qg2);
  QOP_destroy_G(qg);
  if(qg2!=qg) QOP_destroy_G(qg2);
  QOP_asqtad_extract_L_to_qdp(qu, qul, fla);
  for(int i=0; i<4; i++) {
    packM(l, fl[i]->f, qu[i]);
    QDP_destroy_M(qu[i]);
    if(ll) {
      packM(l, ll[i]->f, qul[i]);
      QDP_destroy_M(qul[i]);
    }
  }
}

void
qopDslash(Layout *l, real *x, real *u[8], real *u3[8], real mass, real *y, char *sub)
{
  QDP_ColorMatrix *qu[4], *qu3[4];
  QDP_ColorVector *out, *in;
  in = QDP_create_V();
  out = QDP_create_V();
  unpackV(l, in, y);
  unpackV(l, out, x);
  for(int i=0; i<4; i++) {
    qu[i] = QDP_create_M();
    unpackM(l, qu[i], u[2*i]);
    QLA_Real two = 2;
    QDP_M_eq_r_times_M(qu[i], &two, qu[i], QDP_all);
    if(u3) {
      qu3[i] = QDP_create_M();
      unpackM(l, qu3[i], u3[2*i]);
      QDP_M_eq_r_times_M(qu3[i], &two, qu3[i], QDP_all);
    }
  }
  QOP_FermionLinksAsqtad *fla;
  if(u3) {
    fla = QOP_asqtad_create_L_from_qdp(qu, qu3);
  } else {
    fla = QOP_asqtad_create_L_from_qdp(qu, NULL);
  }
  QOP_evenodd_t eoOut=QOP_EVENODD, eoIn=QOP_EVENODD;
  if(sub[0]=='e') {
    eoOut = QOP_EVEN;
    eoIn = QOP_ODD;
  }
  if(sub[0]=='o') {
    eoOut = QOP_ODD;
    eoIn = QOP_EVEN;
  }
  QOP_asqtad_dslash_qdp(NULL, fla, mass, out, in, eoOut, eoIn);
  QLA_Real n2;
  QDP_r_eq_norm2_V(&n2, out, QDP_all);
  printf0("out2: %g\n", n2);
  packV(l, x, out);
  QDP_destroy_V(in);
  QDP_destroy_V(out);
  for(int i=0; i<4; i++) {
    QDP_destroy_M(qu[i]);
    if(u3) QDP_destroy_M(qu3[i]);
  }
}

void
qopSolveMulti(Layout *l, real *x[], real *u[8], real *u3[8],
	      double masses[], real *y, int nmasses, double rsq, char *sub)
{
  QDP_ColorMatrix *qu[4], *qu3[4];
  QDP_ColorVector *out[nmasses], *in, **outp;
  outp = out;
  in = QDP_create_V();
  unpackV(l, in, y);
  for(int i=0; i<nmasses; i++) {
    out[i] = QDP_create_V();
    unpackV(l, out[i], x[i]);
    QDP_V_eq_zero(out[i], QDP_even);
  }
  for(int i=0; i<4; i++) {
    qu[i] = QDP_create_M();
    unpackM(l, qu[i], u[2*i]);
    QLA_Real two = 2;
    QDP_M_eq_r_times_M(qu[i], &two, qu[i], QDP_all);
    if(u3) {
      qu3[i] = QDP_create_M();
      unpackM(l, qu3[i], u3[2*i]);
      QDP_M_eq_r_times_M(qu3[i], &two, qu3[i], QDP_all);
    }
  }
  QOP_FermionLinksAsqtad *fla;
#ifdef NAIK
  fla = QOP_asqtad_create_L_from_qdp(qu, qu3);
#else
  fla = QOP_asqtad_create_L_from_qdp(qu, NULL);
#endif
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
  QOP_asqtad_invert_multi_qdp(&info, fla, &inv_arg, &rap,
			      &mfp, &nmasses, &outp, &in, 1);
  //QLA_Real n2;
  //QDP_r_eq_norm2_V(&n2, (QDP_ColorVector*)out, QDP_all);
  printf0("QOP its: %i\n", res_arg.final_iter);
  QDP_destroy_V(in);
  for(int i=0; i<nmasses; i++) {
    packV(l, x[i], out[i]);
    QDP_destroy_V(out[i]);
  }
  for(int i=0; i<4; i++) {
    QDP_destroy_M(qu[i]);
    if(u3) QDP_destroy_M(qu3[i]);
  }
}

void
qopSolve(Layout *l, real *x, real *u[8], real *u3[8], real mass, real *y,
	 double rsq, char *sub)
{
  QDP_ColorMatrix *qu[4], *qu3[4];
  QDP_ColorVector *out, *in;
  in = QDP_create_V();
  out = QDP_create_V();
  unpackV(l, in, y);
  unpackV(l, out, x);
  for(int i=0; i<4; i++) {
    qu[i] = QDP_create_M();
    unpackM(l, qu[i], u[2*i]);
    QLA_Real two = 2;
    QDP_M_eq_r_times_M(qu[i], &two, qu[i], QDP_all);
    if(u3) {
      qu3[i] = QDP_create_M();
      unpackM(l, qu3[i], u3[2*i]);
      QDP_M_eq_r_times_M(qu3[i], &two, qu3[i], QDP_all);
    }
  }
  QOP_FermionLinksAsqtad *fla;
#ifdef NAIK
  fla = QOP_asqtad_create_L_from_qdp(qu, qu3);
#else
  fla = QOP_asqtad_create_L_from_qdp(qu, NULL);
#endif
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

  QDP_V_eq_zero(out, QDP_even);
  //QOP_verbose(3);
  QOP_asqtad_invert_qdp(&info, fla, &inv_arg, &res_arg, mass, out, in);
  //QLA_Real n2;
  //QDP_r_eq_norm2_V(&n2, (QDP_ColorVector*)out, QDP_all);
  printf0("QOP its: %i\n", res_arg.final_iter);
  packV(l, x, out);
  QDP_destroy_V(in);
  QDP_destroy_V(out);
  for(int i=0; i<4; i++) {
    QDP_destroy_M(qu[i]);
    if(u3) QDP_destroy_M(qu3[i]);
  }
}

