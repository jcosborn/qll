#include <vec_types.h>
#include <vec_ops.h>

#ifdef DOUBLE
#define P(x) CAT(x,D)
#else
#define P(x) CAT(x,F)
#endif

void P(X_veq_perm_xX)(real *rr, int perm, int *idx, real *vv, int nelem, int off, int len);
void P(xX_veq_yX)(int *idx1, real *rr, int *idx2, real *vv, int nelem, int off, int len);
void P(xX_veq_X)(int *idx1, real *rr, real *vv, int nelem, int off, int len);
void P(X_veq_xXc)(real *rr, int *idx, real *vv, int nelem, int off, int len);
void P(X_veq_xXn)(real *rr, int *idx, real *vv, int nelem, int off, int len);
void P(X_eq_xXp)(real *rr, int *idx1, int perm, real *vv,
		 int nelem, int i0, int i1);
void P(X_veq_xXp)(real *rr, int *idx1, int perm, int *idx2, real *vv,
		  int nelem, int off, int len);
void P(X_veq_xXp2)(real *rr, int *idx1, int perm, real *vv,
		   int nelem, int off, int len);

void P(V_veq_perm_xV)(real *rr, int perm, int *idx, real *vv, int off, int len);
void P(V_veq_r_times_V)(real *rr, real aa, real *vv, int off, int len);
void P(V_vpeq_M_times_xV)(real *rr, real *mm, int *idx, real *vv,
			  int off, int len);
void P(V_vpeq_M_times_xVc)(real *rr, real *mm, int *idx,
			   real *vv, int off, int len);
void P(V_vpeq_M_times_xVp)(real *rr, real *mm, int *idx,
			   int perm, real *vv, int i0, int i1);
void P(V_vpeq_M_times_xVn)(real *rr, real *mm, int *idx,
			   real *vv, int off, int len);
void P(V_vpeq_M_times_xVb)(real *rr, real *mm, int blend, int *idx1, real *vv1,
			   int *idx2, real *vv2, int i0, int i1);
void P(xV_vpeq_M_times_V)(int *idx1, real *rr, real *mm, real *vv, int off, int len);
void P(xV_vpeq_xM_times_V)(int *idx, real *rr, real *mm, real *vv, int off, int len);
void P(xV_vpeq_xM_times_yV)(int *idx1, real *rr, real *mm, int *idx2, real *vv, int off, int len);
void P(xV_vpeq_yM_times_yV)(int *idx1, real *rr, int *idx2, real *mm, real *vv, int begin, int end);
void P(xV_vpeq_xM_times_yVzVb)(int *idx0, real *rr, real *mm, int blend, int *idx1, real *vv1, int *idx2, real *vv2, int i0, int i1);

void P(X_eq_pack_xX)(real *rr, int pack, int *idx, real *vv,
		     int nelem, int i0, int i1);
void P(X_veq_pack_xX)(real *rr, int pack, int *idx, real *vv,
		      int nelem, int off, int len);
void P(xX_veq_blend_xXX)(int *idx1, real *rr, int blend, int *idx2,
			 real *vv1, real *vv2, int nelem, int off, int len);
void P(xX_veq_blend_xXxX)(int *idx0, real *rr, int blend, int *idx1, real *vv1,
			  int *idx2, real *vv2, int nelem, int off, int len);
void P(xX_eq_blend_xXxX)(int *idx0, real *rr, int blend, int *idx1, real *vv1,
			 int *idx2, real *vv2, int nelem, int i0, int i1);
void P(X_veq_blend_xXxXn)(real *rr, int blend, int *idx1, real *vv1,
			  int *idx2, real *vv2, int nelem, int off, int len);

double P(sumElem_X)(real *vv, int nelem, int i0, int i1);
double P(norm2_X)(real *vv, int nelem, int i0, int i1);
void P(r_veq_norm2_X)(double *r, real *vv, int nelem, int i0, int i1);
double P(re_X_dot_X)(real *xx, real *yy, int nelem, int i0, int i1);
void P(r_veq_re_X_dot_X)(double *r, real *vv, real *ww, int nelem, int i0, int i1);
void P(X_veq_r)(real *rr, real aa, int nelem, int i0, int i1);
void P(X_veq_r_times_X)(real *rr, real aa, real *vv,
			int nelem, int i0, int i1);
void P(X_vpeq_r_times_X)(real *rr, real aa, real *vv,
			 int nelem, int i0, int i1);
void P(X_veq_X_plus_r_times_X)(real *rr, real *vv, real aa, real *ww,
			       int nelem, int i0, int i1);
void P(X_eq_scale_plus_X)(real *rr, real aa, real *vv, int nelem, int i0, int i1);
void P(X_eq_scale_plus_r_times_X)(real *rr, real aa, real bb, real *vv, int nelem, int i0, int i1);
void P(V_vmeq_xMap_times_xVp)(real *rr, real *mm, int *idx,
			      int perm, real *vv, int i0, int i1);
void P(X_peq_M_times_X)(real *rr, real *mm, real *vv, int nelem, int nc,
			int tdv, int i0, int i1);
void P(X_peq_M_times_Xc)(real *rr, real *mm, real *vv, int nelem, int nc,
			 int tdv, int i0, int i1);
void P(X_peq_M_times_Xa)(real *rr, real *mm, real *vv, int nelem, int nc,
			 int tdv, int i0, int i1);
void P(X_peq_Ma_times_X)(real *rr, real *mm, real *vv, int nelem, int nc,
			 int tdv, int i0, int i1);
void P(X_peq_M_times_xXp)(real *rr, real *mm, int *idx, int perm, real *vv,
			  int nelem, int nc, int tdv, int i0, int i1);
void P(xX_peq_xM_times_yXzXb)(int *idx0, real *rr, real *mm, int blend,
			      int *idx1, real *vv1, int *idx2, real *vv2,
			      int nelem, int nc, int tdv, int i0, int i1);
void P(X_eq_Xa)(real *rr, real *vv, int nr, int nc, int i0, int i1);
void P(Xadj)(real *rr, int nc, int i0, int i1);

void P(makeU)(real *mm, int nc, int i0, int i1);

void
P(H_eq_spinproj_D)(real *rr, real *vv, int dir, int sign, int i0, int i1);
void
P(D_eq_spinrecon_H)(real *rr, real *vv, int dir, int sign, int i0, int i1);

void
P(D_peq_spinproj_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
			      real *vv, int dir, int sign, int i0, int i1);
void
P(D_meq_spinproj_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
			      real *vv, int dir, int sign, int i0, int i1);

void P(D_meq_spinproj0p_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
				     real *vv, int i0, int i1);
void P(D_meq_spinproj0m_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
				     real *vv, int i0, int i1);
void P(D_meq_spinproj1p_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
				     real *vv, int i0, int i1);
void P(D_meq_spinproj1m_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
				     real *vv, int i0, int i1);
void P(D_meq_spinproj2p_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
				     real *vv, int i0, int i1);
void P(D_meq_spinproj2m_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
				     real *vv, int i0, int i1);
void P(D_meq_spinproj3p_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
				     real *vv, int i0, int i1);
void P(D_meq_spinproj3m_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
				     real *vv, int i0, int i1);
void P(D_meq_spinproj4p_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
				     real *vv, int i0, int i1);
void P(D_meq_spinproj4m_M_times_xDp)(real *rr, real *mm, int *idx, int perm,
				     real *vv, int i0, int i1);

#define X_eq_r X_veq_r
#define X_veq_r P(X_veq_r)
#define X_veq_r_times_X P(X_veq_r_times_X)
#define X_vpeq_r_times_X P(X_vpeq_r_times_X)
#define r_veq_norm2_X P(r_veq_norm2_X)
#define sumElem_X P(sumElem_X)
#define norm2_X P(norm2_X)
#define re_X_dot_X P(re_X_dot_X)
#define r_veq_re_X_dot_X P(r_veq_re_X_dot_X)
#define X_veq_X_plus_r_times_X P(X_veq_X_plus_r_times_X)
#define X_eq_scale_plus_X P(X_eq_scale_plus_X)
#define X_eq_scale_plus_r_times_X P(X_eq_scale_plus_r_times_X)
#define X_eq_Xa P(X_eq_Xa)
#define Xadj P(Xadj)

#define V_eq_zero(x) X_veq_r(x, 0, 1, ti0, ti1)
#define V_eq_V(x,y)  X_veq_r_times_X(x, 1, y, 1, ti0, ti1)
#define r_eq_norm2_V(r,x) r_veq_norm2_X(r,x,1,ti0,ti1)
#define r_eq_re_V_dot_V(r,x,y) r_veq_re_X_dot_X(r,x,y,1,ti0,ti1)
#define V_peq_r_times_V(x,a,y) X_veq_X_plus_r_times_X(x,x,a,y,1,ti0,ti1)
#define V_meq_r_times_V(x,a,y) X_veq_X_plus_r_times_X(x,x,-(a),y,1,ti0,ti1)
#define V_eq_scale_plus_V(x,a,y) X_eq_scale_plus_X(x,a,y,1,ti0,ti1)
#define V_eq_scale_plus_r_times_V(x,a,b,y) X_eq_scale_plus_r_times_X(x,a,b,y,1,ti0,ti1)
