#-----------------------------------------------------------------------
#
#                     Generic GSLIB Makefile
#                     **********************
#
progname=sisim
#
#
#-----------------------------------------------------------------------


PROG   = ../bin/$(progname)
FFLAGS = -g
SRCS   =  $(progname).f
OBJS   =  $(progname).o
INCS   =  $(progname).inc
LIBS   = ../gslib/libgs.a
all:$(PROG) 
$(PROG):$(OBJS)
	f90 -o $@ $(OBJS) $(LIBS)
.f.o:
	f90 -c $(FFLAGS) $<
$(OBJS):$(INCS)
clean:
	/bin/rm -rf $(OBJS) $(PROG)
