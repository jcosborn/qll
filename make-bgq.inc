#TAG = -icc
#TAG = -icc-omp
#TAG = -mic
SCIDAC = $(HOME)/lqcd/install
QOPDIR = $(SCIDAC)/qopqdp-spi$(TAG)
QDPDIR = $(SCIDAC)/qdp-spi$(TAG)
QLADIR = $(SCIDAC)/qla$(TAG)
QIODIR = $(SCIDAC)/qio$(TAG)
QMPDIR = $(SCIDAC)/qmp-spi$(TAG)
#QMPDIR = $(SCIDAC)/qmp$(TAG)
QMPCONF = $(QMPDIR)/bin/qmp-config
QLACONF = $(QLADIR)/bin/qla-config

#OMP = -fopenmp
OMP = -qsmp=omp

COPT = $(shell $(QLACONF) --copts) -g -DQDP_PROFILE
#COPT = -O3 -fopenmp -g -DQDP_PROFILE
#COPT += -DC1
#COPT += -DC2
#COPT += -DC4
COPT += -DQPX -DTIBDEFAULT=4
#COPT += -DAVX

CC = $(shell $(QMPCONF) --cc)
#CC = tmpwrap /soft/perftools/autoperf/scripts/xl.legacy.ndebug/mpixlc_r
#CC = tmpwrap mpicc
#LD = tmpwrap hpclink mpixlc_r
