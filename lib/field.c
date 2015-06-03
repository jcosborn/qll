#include <common.h>

void
fieldSet(Field *f, real *r, Layout *l, int nelem/*, int nc*/)
{
  f->l = l;
  f->f = r;
  f->nelem = nelem;
  //f->nc = nc;
}

void
fieldNew(Field *f, Layout *l, int nelem/*, int nc*/)
{
  f->l = l;
  f->f = myalloc(nelem*l->nSites*sizeof(real));
  f->nelem = nelem;
  //f->nc = nc;
}

void
fieldFree(Field *f)
{
  free(f->f);
}

void
fieldZero(Field *f, Subset *sub)
{
  BEGIN_PARALLEL
    {
      SUBSET_DEFS2;
      setOffSub2(*sub);
      X_eq_r(f->f, 0, f->nelem, ti0, ti1);
    }
  END_PARALLEL;
}

double
fieldSumElem(Field *f, Subset *sub)
{
  double n2 = 0;
  beginParallelSubset(*sub) {
    double r;
    r = sumElem_X(f->f, f->nelem, ti0, ti1);
    sumDoubleArray(&r, 1);
    singleThread { n2 = r; }
  } endParallelSubset;
  return n2;
}

double
fieldNorm2(Field *f, Subset *sub)
{
  double n2 = 0;
  beginParallelSubset(*sub) {
    double r;
    r = lnorm2_F(f);
    sumDoubleArray(&r, 1);
    singleThread { n2 = r; }
  } endParallelSubset;
  return n2;
}

#if 0
void
fieldLocalNorm2(real *r, Field *f, Subset *sub)
{
  BEGIN_PARALLEL
    {
      SUBSET_DEFS2;
      setOffSub2(*sub);
      X_eq_r(f->f, 0, f->nelem, ti0, ti1);
    }
  END_PARALLEL;
}
#endif
