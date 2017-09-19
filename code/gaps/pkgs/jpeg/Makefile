#
# Package name
#

NAME=jpeg



#
# List of source files
#

# library object files common to compression and decompression
COMSRCS= jcomapi.c jutils.c jerror.c jmemmgr.c jmemnobs.c

# compression library object files
CLIBSRCS= jcapimin.c jcapistd.c jctrans.c jcparam.c \
        jdatadst.c jcinit.c jcmaster.c jcmarker.c jcmainct.c \
        jcprepct.c jccoefct.c jccolor.c jcsample.c jchuff.c \
        jcphuff.c jcdctmgr.c jfdctfst.c jfdctflt.c \
        jfdctint.c

# decompression library object files
DLIBSRCS= jdapimin.c jdapistd.c jdtrans.c jdatasrc.c \
        jdmaster.c jdinput.c jdmarker.c jdhuff.c jdphuff.c \
        jdmainct.c jdcoefct.c jdpostct.c jddctmgr.c \
        jidctfst.c jidctflt.c jidctint.c jidctred.c \
        jdsample.c jdcolor.c jquant1.c jquant2.c jdmerge.c

# JPEG library source files
CSRCS=$(CLIBSRCS) $(DLIBSRCS) $(COMSRCS)



#
# GAPS library makefile
#

include ../../makefiles/Makefile.pkgs



