##
##	PDP1 Makefile 
##
FILE_EXT = linux
##
XGPATH =  /home/aerosun/a/share/lib
##	Gas type for the mcc package.  Uncomment only one.
##
##
##
EXEC = xpdp1.$(FILE_EXT)
##
CC = gcc
##
##      C compiler.  Normally gcc.
CFLAGS= -g -I$(XGPATH) -DPDP1PATH=\"$(HOME)/Birdsall/xpdp1/src\" 
##
##      Flags used for compiling.  -O2 for optimization, which is optional.
##	If the X11 include files are not located in /usr/include/X11, you
##	must specify the directory as a flag prefixed by -I.  For instance, 
##	if the include files for X11 are located in /usr/local/include, the 
##	flag -I/usr/local/include must also be passed.
##
#LIBS  = -L$(XGPATH) -L/usr/local/lib -L/usr/X11R6/lib -lXGC250 -ltk -ltcl -lXpm -lX11 -lm -ldl 
##
MPI_INCLUDE = -I/opt/mpich2-intel/1.4.1p1/include
MPI_LIBDIR =-L/opt/mpich2-intel/1.4.1p1/lib64
MPI_LIB =-lmpich

X11_LIBDIR =-L/usr/X11R6/lib64/
X11_INCLUDE =-I/usr/X11R6/lib64/include
X11_LIB = -lX11

XGRAFIX_LIBDIR =-L/home/aerosun/a/share/lib
XGRAFIX_INCLUDE =-I/home/aerosun/a/share/include
XGRAFIX_LIB = -lXGC250

TCL_LIBDIR = -L/lib64
TCL_LIB = -ltk -ltcl
XPM_LIB = -lXpm

LIB64_DIR = -L/usr/lib64
LOCAL_LIBDIR = -L/usr/local/lib

LIBDIRS = $(MPI_LIBDIR) $(X11_LIBDIR) $(XGRAFIX_LIBDIR) $(TCL_LIBDIR) $(LIB64_DIR) $(LOCAL_LIBDIR)

LIBS = $(MPI_LIB) $(X11_LIB )$(XGRAFIX_LIB) $(TCL_LIB) -lm -ldl -L/home/aerosun/a/share/lib -L/usr/lib64 -L/home/aerosun/a/vayyaswa/

INCLUDES = $(MPI_INCLUDE) $(X11_INCLUDE) $(XGRAFIX_INCLUDE)

##
## 	Libraries and their directories used in loading.  Normal is -lm and -lX11
##	for the Math and X11 libraries, which is NOT optional.  On Unicos, -lnet
##	is also used.  If the X11 libraries are not located in /usr/lib/X11, you
##	must also specify the directory as a flag prefixed by -L.  For instance,
##	if the include files for X11 are located in /usr/local/lib, the flag
##	-L/usr/local/lib must also be passed.  On Sun, use -Bstatic to avoid
##	problems with differing versions of OS.  Man cc for more information on flags.
##
##
PDP1OBJ= fft.o field.o move.o gather.o pdp1.o start.o load.o prest.o \
	 padjus.o initwin.o maxwellv.o  mccdiaginit.o xsect.o\
	 argonmcc.o heliummcc.o neonmcc.o oxygenmcc.o nitrogen.o mcc.o

all:	$(PDP1OBJ) $(EXEC)

.c.o:	pdp1.h xsect.h
	$(CC) -c $(CFLAGS) $(INCLUDES) $*.c

$(EXEC):	$(PDP1OBJ)
		$(CC) -o $(EXEC) $(PDP1OBJ) $(LIBDIRS) $(LIBS)

clean:
	rm *.o *~
