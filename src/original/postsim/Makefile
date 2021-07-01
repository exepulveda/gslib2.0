#-----------------------------------------------------------------------
#
#                     Generic GSLIB Makefile
#                     **********************
#
progname=postsim
#
#
#-----------------------------------------------------------------------


PROG   = ../bin/$(progname)
FFLAGS = -O
SRCS   =  $(progname).f
OBJS   =  $(progname).o
INCS   =  
LIBS   = ../gslib/libgs.a
all:$(PROG) 
$(PROG):$(OBJS)
	f90 -o $@ $(OBJS) $(LIBS)
.f.o:
	f90 -c $(FFLAGS) $<
$(OBJS):$(INCS)
clean:
	/bin/rm -rf $(OBJS) $(PROG)
