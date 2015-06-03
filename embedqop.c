#include <qop.h>
#include <qop_internal.h>
#include <qll.h>

#if QOP_Precision == 'F'
#define P(x) x ## F
typedef float real;
#else
#define P(x) x ## D
typedef double real;
#endif

void
get_rankGeom(QDP_Lattice *lat, int myrank, int nd, int *ls, int *rg)
{
  int n = QDP_numsites_L(lat, myrank);
  int xmin[nd], xmax[nd];
  for(int i=0; i<nd; i++) {
    xmin[i] = ls[i];
    xmax[i] = 0;
  }
  for(int i=0; i<n; i++) {
    int x[nd];
    QDP_get_coords(x, myrank, i);
    for(int j=0; j<nd; j++) {
      if(x[j]<xmin[j]) xmin[j] = x[j];
      if(x[j]>xmax[j]) xmax[j] = x[j];
    }
  }
  for(int i=0; i<nd; i++) {
    rg[i] = ls[i]/(1+xmax[i]-xmin[i]);
  }
}

static int setup = 0;
static Layout layout;
static P(StaggeredSolver) sse;
static QLA_Real **u=NULL, **u3=NULL;
static int nl = 0;

void
setup_qll(QDP_Lattice *lat)
{
  if(!setup) {
    setup = 1;
    layout.nranks = QMP_get_number_of_nodes();
    layout.myrank = QMP_get_node_number();
    int nd = QDP_ndim_L(lat);
    int ls[nd], rg[nd];
    QDP_latsize_L(lat, ls);
    get_rankGeom(lat, layout.myrank, nd, ls, rg);
    P(stagDslashSetup)(&layout, nd, ls, rg);
  }
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
#pragma omp parallel for
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
#pragma omp parallel for
  for(int j=0; j<nsites; j++) { // qll sites
    int x[nd];
    LayoutIndex li;
    li.rank = l->myrank;
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
setup_qll_solver(QOP_FermionLinksAsqtad *fla)
{
  if(!fla->dblstored) {
    if(layout.myrank==0) printf("error: links not double stored\n");
  }
  nl = fla->nlinks;
  //printf("nl: %i\n", nl);
  if(nl==8) { // no naik
    u = myalloc(nl*sizeof(QLA_Real*));
    u3 = NULL;
    for(int i=0; i<nl; i++) {
      u[i] = myalloc(layout.nSites*18*sizeof(QLA_Real));
      packM(&layout, u[i], fla->dbllinks[i]);
    }
  } else { // naik
    nl /= 2;
    u = myalloc(nl*sizeof(QLA_Real*));
    u3 = myalloc(nl*sizeof(QLA_Real*));
    for(int i=0; i<nl; i++) {
      u[i] = myalloc(layout.nSites*18*sizeof(QLA_Real));
      u3[i] = myalloc(layout.nSites*18*sizeof(QLA_Real));
      packM(&layout, u[i], fla->dbllinks[2*i]);
      packM(&layout, u3[i], fla->dbllinks[2*i+1]);
    }
  }
  P(createStaggeredSolver)(&sse, &layout, u, u3, "even");
  //sse.noscale = 1;
}

void
free_qll_solver(void)
{
  P(freeStaggeredSolver)(&sse);
  if(u) {
    for(int i=0; i<nl; i++) {
      free(u[i]);
    }
    free(u);
    u = NULL;
  }
  if(u3) {
    for(int i=0; i<nl; i++) {
      free(u3[i]);
    }
    free(u3);
    u3 = NULL;
  }
}

void
solve_qll(QDP_ColorVector *dest, QDP_ColorVector *src, double mass,
	  QOP_invert_arg_t *inv_arg, QOP_resid_arg_t *res_arg)
{
  QLA_Real *s = myalloc(layout.nSites*6*sizeof(QLA_Real));
  QLA_Real *d = myalloc(layout.nSites*6*sizeof(QLA_Real));
  packV(&layout, s, src);
  packV(&layout, d, dest);
  double rsq = res_arg->rsqmin;
  int maxits = inv_arg->max_iter;
  int its = P(solve2)(&sse, d, mass, s, rsq, maxits, "even");
  res_arg->final_iter = its;
  unpackV(&layout, dest, d);
  free(s);
  free(d);
}

int
main(int argc, char *argv[])
{
  QDP_initialize(&argc, &argv);
  QDP_verbose(0);
  QDP_profcontrol(0);

  int nd = 4;
  if(argc>1) nd = argc - 1;
  int ls[nd];
  if(argc>1) {
    for(int i=0; i<nd; i++) ls[i] = atoi(argv[i+1]);
  } else {
    for(int i=0; i<nd; i++) ls[i] = 4;
  }
  QDP_set_latsize(nd, ls);
  QDP_create_layout();
  QDP_ColorMatrix *u[nd];
  for(int i=0; i<nd; i++) {
    QLA_Complex z;
    QLA_c_eq_r(z, 1);
    u[i] = QDP_create_M();
    QDP_M_eq_c(u[i], &z, QDP_all);
  }
  QDP_ColorVector *src, *dest1, *dest2;
  src = QDP_create_V();
  dest1 = QDP_create_V();
  dest2 = QDP_create_V();
  QDP_V_eq_zero(src, QDP_all);
  QDP_V_eq_zero(dest1, QDP_all);
  QDP_V_eq_zero(dest2, QDP_all);
  if(QDP_this_node==0) {
    QLA_Real *x = QDP_site_ptr_readwrite_V(src, 0);
    *x = 1;
  }

  QOP_layout_t qoplayout = QOP_LAYOUT_ZERO;
  qoplayout.latdim = nd;
  qoplayout.latsize = (int *) myalloc(nd*sizeof(int));
  for(int i=0; i<nd; i++) {
    qoplayout.latsize[i] = ls[i];
  }
  qoplayout.machdim = -1;
  QOP_init(&qoplayout);

  QOP_info_t info = QOP_INFO_ZERO;
  QOP_asqtad_coeffs_t coeffs = QOP_ASQTAD_COEFFS_ZERO;
  coeffs.one_link = 1;
  coeffs.three_staple = 0.1;
  coeffs.five_staple = 0.1;
  coeffs.seven_staple = 0.1;
  coeffs.lepage = 0.1;
  coeffs.naik = 0.1;
  QOP_invert_arg_t inv_arg = QOP_INVERT_ARG_DEFAULT;
  inv_arg.max_iter = 600;
  inv_arg.restart = 200;
  inv_arg.max_restarts = 5;
  inv_arg.evenodd = QOP_EVEN;
  inv_arg.mixed_rsq = 0;
  QOP_resid_arg_t res_arg = QOP_RESID_ARG_DEFAULT;
  res_arg.rsqmin = 1e-4;

  double mass = 0.9;

  QOP_opt_t optcg;
  optcg.tag = "cg";
  optcg.value = 1;
  QOP_asqtad_invert_set_opts(&optcg, 1);
  QOP_GaugeField *gf = QOP_convert_G_from_qdp(u);
  QOP_FermionLinksAsqtad *fla = QOP_asqtad_create_L_from_G(&info, &coeffs, gf);

  QLA_Real nrm2;
  QDP_r_eq_norm2_V(&nrm2, src, QDP_all);
  printf("src2: %g\n", nrm2);
  QOP_asqtad_invert_qdp(&info, fla, &inv_arg, &res_arg, mass, dest1, src);
  QDP_r_eq_norm2_V(&nrm2, dest1, QDP_all);
  printf("dest12: %g\n", nrm2);
  
  QDP_Lattice *lat = QDP_get_lattice_V(src);
  //for(int i=0; i<1000; i++) {
    setup_qll(lat);
    setup_qll_solver(fla);
    solve_qll(dest2, src, mass, &inv_arg, &res_arg);
    free_qll_solver();
    //sleep(1);
    //}

  QDP_V_meq_V(dest1, dest2, QDP_all);
  QDP_r_eq_norm2_V(&nrm2, dest1, QDP_all);
  printf("diff2: %g\n", nrm2);

  QDP_finalize();
}
