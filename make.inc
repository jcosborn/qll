#TAG = -icc
#TAG = -icc-omp
#TAG = -mic
SCIDAC = $(HOME)/lqcd/install
QOPDIR = $(SCIDAC)/qopqdp$(TAG)
QDPDIR = $(SCIDAC)/qdp$(TAG)
QLADIR = $(SCIDAC)/qla$(TAG)
QIODIR = $(SCIDAC)/qio$(TAG)
QMPDIR = $(SCIDAC)/qmp$(TAG)
QMPCONF = $(QMPDIR)/bin/qmp-config
QLACONF = $(QLADIR)/bin/qla-config

OMP = -fopenmp
#OMP = -qsmp=omp

#OPT2 = -opt-prefetch-distance=64,8 -opt-streaming-cache-evict=0 -opt-streaming-stores always

#COPT = $(shell $(QLACONF) --copts) -g -DQDP_PROFILE $(OPT2)
COPT = -Wall -std=c99 -O3 -march=core-avx2 -ffast-math -funroll-loops -fprefetch-loop-arrays -fpeel-loops -funswitch-loops -g -DQDP_PROFILE
#COPT += -DC1
#COPT += -DC2
#COPT += -DC4
#COPT += -DC8
COPT += -DAVX -DUSE_GCC_VECTOR_OPS
#COPT += -DMIC

CC = $(shell $(QMPCONF) --cc)
