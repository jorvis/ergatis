BUILD=Release

###  PATH TO A PRE-BUILT C++ TOOLKIT
include Makefile.builddir


###  DEFAULT COMPILATION FLAGS  -- DON'T EDIT OR MOVE THESE 4 LINES !!!
include $(builddir)/Makefile.mk
srcdir = .
BINCOPY = @:
LOCAL_CPPFLAGS = -I.


#############################################################################
###  EDIT SETTINGS FOR THE DEFAULT (APPLICATION) TARGET HERE              ### 
APP = Bsml2Asn1 
SRC = Bsml2Asn1 Bsml2Asn1Converter FastaIndexer Registry
# OBJ =

# PRE_LIBS = $(NCBI_C_LIBPATH) .....

LIB = submit taxon1 xobjutil xobjmgr $(OBJMGR_LIBS)
LIBS = $(CMPRS_LIBS) $(NETWORK_LIBS) $(DL_LIBS) $(ORIG_LIBS)


## If you need the C toolkit...
# LIBS     = $(NCBI_C_LIBPATH) $(NCBI_C_ncbi) $(NETWORK_LIBS) $(ORIG_LIBS)
# CPPFLAGS = $(ORIG_CPPFLAGS) $(NCBI_C_INCLUDE)

# CFLAGS   = $(ORIG_CFLAGS)
# CXXFLAGS = $(ORIG_CXXFLAGS)
# LDFLAGS  = $(ORIG_LDFLAGS)
#                                                                         ###
#############################################################################


###  APPLICATION BUILD RULES  -- DON'T EDIT OR MOVE THIS LINE !!!
include $(builddir)/Makefile.app
MAKEFILE = Makefile.test_app


###  PUT YOUR OWN ADDITIONAL TARGETS (MAKE COMMANDS/RULES) HERE
XERCES = /usr/local/devel/ANNOTATION/xerces/xerces-c-src_2_7_0
XERCES_LIB_PATH = $(XERCES)/lib
XERCES_INC_PATH = $(XERCES)/include
XALAN = /usr/local/devel/ANNOTATION/xalan/xml-xalan/c
XALAN_LIB_PATH = $(XALAN)/lib
XALAN_INC_PATH= $(XALAN)/src

LDFLAGS = $(ORIG_LDFLAGS) -Wl,-rpath,$(XERCES_LIB_PATH) -lxerces-c -L$(XERCES_LIB_PATH) -Wl,-rpath,$(XALAN_LIB_PATH) -lxalan-c -lxalanMsg -L$(XALAN_LIB_PATH)

ifeq ($(BUILD), Debug)
TMP_CPPFLAGS = $(shell echo '$(ORIG_CPPFLAGS)' | sed -e 's/-DNDEBUG/-D_DEBUG/')
TMP_CXXFLAGS = $(shell echo '$(ORIG_CFLAGS)' | sed -e 's/-O./-g/')
else
TMP_CPPFLAGS = $(ORIG_CPPFLAGS)
TMP_CXXFLAGS = $(ORIG_CXXFLAGS)
endif

CPPFLAGS = $(TMP_CPPFLAGS) -I$(XERCES_INC_PATH) -I$(XALAN_INC_PATH)
CXXFLAGS = $(TMP_CXXFLAGS)

#CPPFLAGS = $(ORIG_CPPFLAGS) -I$(XERCES_INC_PATH) -I$(XALAN_INC_PATH)
