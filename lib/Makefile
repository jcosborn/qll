include ../make.inc

CFLAGS = -I../include

IQMP = $(shell $(QMPCONF) --cflags)
LQMP = $(shell $(QMPCONF) --ldflags) $(shell $(QMPCONF) --libs)
CFLAGSQ = $(COPT) -I../include $(IQMP)
CFLAGS = $(CFLAGSQ) $(OMP)

OBJS = \
common.o \
gather.o \
layout.o \
shifts.o \
vlaD.o \
vlaF.o \
vlawD.o \
vlawF.o \
fieldD.o \
fieldF.o \
gaugeD.o \
gaugeF.o \
transporterF.o \
transporterD.o \
stagDslashF.o \
stagDslashD.o \
wilsonDslashF.o \
wilsonDslashD.o \
cgF.o \
cgD.o \
cgmsF.o \
cgmsD.o

LIB = libqll.a

lib: $(LIB)
$(LIB): $(OBJS)
	ar cru $@ $^

vla0.c: vla0.c.lua
	lua $^ >$@

vlaw0.c: vlaw0.c.lua
	lua $^ >$@

vlaF.o: vlaF.c vla.c vla0.c
vlaD.o: vlaD.c vla.c vla0.c
vlawF.o: vlawF.c vlaw.c vlaw0.c
vlawD.o: vlawD.c vlaw.c vlaw0.c

fieldF.o: fieldF.c field.c
fieldD.o: fieldD.c field.c
gaugeF.o: gaugeF.c gauge.c
gaugeD.o: gaugeD.c gauge.c
transporterF.o: transporterF.c transporter.c
transporterD.o: transporterD.c transporter.c

stagDslashF.o: stagDslashF.c stagDslash.c
stagDslashD.o: stagDslashD.c stagDslash.c
wilsonDslashF.o: wilsonDslashF.c wilsonDslash.c
wilsonDslashD.o: wilsonDslashD.c wilsonDslash.c

cgF.o: cgF.c cg.c
cgD.o: cgD.c cg.c
cgmsF.o: cgmsF.c cgms.c
cgmsD.o: cgmsD.c cgms.c

clean:
	\rm -f $(OBJS) $(LIB)
