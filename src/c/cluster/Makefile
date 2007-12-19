# $Author: bhaas $ 
# $Header: /usr/local/src/cvsroot/ANNOTATION/CLUSTER/src/Makefile,v 1.1 2003/07/07 16:54:09 bhaas Exp $ $Date: 2003/07/07 16:54:09 $ 
#
# $Log: Makefile,v $
# Revision 1.1  2003/07/07 16:54:09  bhaas
# *** empty log message ***
#
# Revision 1.2  1998/07/30 22:50:12  ania
# Rewritten, works fine on Sun.
#
# Revision 1.1  1998/07/30 22:47:57  ania
# Initial revision
#
#
#---------------------------------------------------------------------
#                cluster Makefile 
#
# This makefile can create:
# cluster - program for preparing display coordinates for 
#	   multi-sequence alignments
#
#------------------------------------------------------------------------
SRC_DIR = .
OBJ_DIR = .
BIN_DIR = .

CC = gcc
CFLAGS = -c -I$(OBJ_DIR)  -g
LIBS =  -lm -lnsl

OBJECTS = \
	  $(OBJ_DIR)/cluster.o  \
	  $(OBJ_DIR)/cluster_lists.o  \
	  $(OBJ_DIR)/cluster_inp.o \
	  $(OBJ_DIR)/err.o
SOURCES = \
	$(SRC_DIR)/cluster.c \
	$(SRC_DIR)/cluster_inp.c \
	$(SRC_DIR)/cluster_lists.c \
	$(SRC_DIR)/err.c
INCLUDES = cluster_defs.h \
	  cluster_inp_defs.h \
	  cluster_inp_vars.h \
	  cluster_lists_defs.h \
	  cluster_vars.h

cluster: $(OBJECTS) $(INCLUDES)
	$(CC) $(OBJECTS) $(LIBS) -o $(BIN_DIR)/cluster

$(OBJ_DIR)/cluster.o: $(SRC_DIR)/cluster.c
	$(CC) $(CFLAGS) -c $(SRC_DIR)/cluster.c -o $@

$(OBJ_DIR)/cluster_lists.o: $(SRC_DIR)/cluster_lists.c
	$(CC) $(CFLAGS) -c $(SRC_DIR)/cluster_lists.c -o $@

$(OBJ_DIR)/cluster_inp.o: $(SRC_DIR)/cluster_inp.c
	$(CC) $(CFLAGS) -c $(SRC_DIR)/cluster_inp.c -o $@

$(OBJ_DIR)/err.o: $(SRC_DIR)/err.c
	$(CC) $(CFLAGS) -c $(SRC_DIR)/err.c -o $@

rcs_all:
	co -l cluster.c
	co -l cluster.input
	co -l cluster_defs.h
	co -l cluster_inp.c
	co -l cluster_inp_defs.h
	co -l cluster_inp_vars.h
	co -l cluster_lists.c
	co -l cluster_lists_defs.h
	co -l cluster_vars.h
	co -l err.c

clean:
	rm -f ${OBJ_DIR}/*.o

