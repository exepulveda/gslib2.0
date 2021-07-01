#-----------------------------------------------------------------------
#
#                     Generic GSLIB Makefile
#                     **********************
#
progname=getdata
#
#
#-----------------------------------------------------------------------


PROG   =  $(progname)
FFLAGS = -O
SRCS   =  $(progname).f
OBJS   =  $(progname).o
INCS   =  
LIBS   = 
all:$(PROG) 
$(PROG):$(OBJS)
	f90 -o $@ $(OBJS) $(LIBS)
.f.o:
	f90 -c $(FFLAGS) $<
$(OBJS):$(INCS)
clean:
	/bin/rm -rf $(OBJS) $(PROG)
