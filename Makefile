include make.inc

IQOP = -I$(QOPDIR)/include
LQOP = -L$(QOPDIR)/lib -lqopqdp

IQDP = -I$(QDPDIR)/include
LQDP = -L$(QDPDIR)/lib -lqdp_f -lqdp_d -lqdp_df -lqdp_f3 -lqdp_d3 -lqdp_df3 -lqdp_fn -lqdp_dn -lqdp_dfn -lqdp_int -lqdp_common

IQLA = -I$(QLADIR)/include
LQLA = -L$(QLADIR)/lib -lqla_f -lqla_d -lqla_q -lqla_df -lqla_dq -lqla_f3 -lqla_d3 -lqla_q3 -lqla_df3 -lqla_dq3 -lqla_fn -lqla_dn -lqla_qn -lqla_dfn -lqla_dqn -lqla_int -lqla_random -lqla_cmath

IQIO = -I$(QIODIR)/include
LQIO = -L$(QIODIR)/lib -lqio -llime

IQMP = $(shell $(QMPCONF) --cflags)
LQMP = $(shell $(QMPCONF) --ldflags) $(shell $(QMPCONF) --libs)

CFLAGSQ = $(COPT) -Iinclude $(IQOP) $(IQDP) $(IQLA) $(IQIO) $(IQMP)
CFLAGS = $(CFLAGSQ) $(OMP)
CFLAGSF = $(CFLAGS) -DQOP_PrecisionInt=1 -DQDP_PrecisionInt=1
CFLAGSD = $(CFLAGS) -DQOP_PrecisionInt=2 -DQDP_PrecisionInt=2 -DDOUBLE
LDLIBS = $(LQCD) $(LQOP) $(LQDP) $(LQLA) $(LQIO) $(LQMP) $(LAPACK) -lm $(OMP)

PROGS = testwilsonF testdslashF testdslashD embedqopF embedqopD fieldtestF \
 stressD
OBJSF = qopstagF.o testdslashF.o
OBJSD = qopstagD.o testdslashD.o
WOBJSF = qopwilsonF.o testwilsonF.o
OBJS = $(OBJSF) $(OBJSD) $(WOBJSF) embedqopF.o embedqopD.o fieldtestF.o \
 stressD.o
LIBS = lib/libqll.a
LDFLAGS = $(COPT)

all: $(PROGS)

SUBDIRS = include lib

subdirs: $(SUBDIRS)
lib: include
$(LIBS): lib

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: subdirs $(SUBDIRS)

$(OBJS): include/vla.h include/vec_types.h include/vec_ops.h $(LIBS)

testdslashF: $(OBJSF) $(LIBS)
testdslashD: $(OBJSD) $(LIBS)

embedqopF: embedqopF.o $(LIBS)
embedqopD: embedqopD.o $(LIBS)

testwilsonF: $(WOBJSF) $(LIBS)

fieldtestF: fieldtestF.o qopstagF.o $(LIBS)

stressD: stressD.o $(LIBS)

%F.o: %.c
	$(CC) $(CFLAGSF) -c $< -o $@

%D.o: %.c
	$(CC) $(CFLAGSD) -c $< -o $@

clean: subdirs
	-rm -rf $(PROGS) $(OBJS)
