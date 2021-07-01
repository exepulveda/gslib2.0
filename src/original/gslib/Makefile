#
# Makefile for generating "libgs.a"
#
LIBGS = libgs.a
SRCS = acorni.f backtr.f beyond.f blue.f chknam.f chktitle.f cova3.f \
	dlocate.f dpowint.f dsortem.f gauinv.f gcum.f getindx.f getz.f \
	green.f hexa.f ksol.f ktsol.f locate.f nscore.f numtext.f \
	ordrel.f picksupr.f powint.f psfill.f psline.f pstext.f rand.f \
	red.f resc.f scal.f setrot.f setsupr.f sortem.f sqdist.f \
	srchsupr.f strlen.f
OBJS = $(SRCS:.f=.o)

FFLAGS = -c -O

all: $(LIBGS)

$(LIBGS): $(OBJS)
	ar q $@ $(OBJS)
	ar ts $@

.f.o:
	f90 $(FFLAGS) $*.f

install:
	install -s -m 644 $(LIBGS) $(LIBDIR)

clean:
	/bin/rm -rf $(OBJS) $(LIBGS) core
