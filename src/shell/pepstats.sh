#!/bin/sh
if [ -z "$SCHEMA_DOCS_DIR" ]
then
    SCHEMA_DOCS_DIR=/usr/local/devel/ANNOTATION/ard/testing_automated/docs
export SCHEMA_DOCS_DIR
fi
if [ -z "$WORKFLOW_WRAPPERS_DIR" ]
then
    WORKFLOW_WRAPPERS_DIR=/usr/local/devel/ANNOTATION/ard/testing_automated/bin
export WORKFLOW_WRAPPERS_DIR
fi
if [ -z "$WORKFLOW_DOCS_DIR" ]
then
    WORKFLOW_DOCS_DIR=/usr/local/devel/ANNOTATION/ard/testing_automated/docs
export WORKFLOW_DOCS_DIR
fi

umask 0000

unset PERL5LIB
unset LD_LIBRARY_PATH

LANG=C
export LANG
LC_ALL=C
export LC_ALL

EMBOSS_PATH=/usr/local/packages/EMBOSS
export EMBOSS_PATH

PATH=$EMBOSS_PATH/bin:$PATH
export PATH
EMBOSS_ACDROOT=$EMBOSS_PATH/share/EMBOSS/acd
export EMBOSS_ACDROOT
EMBOSS_DATA=$EMBOSS_PATH/share/EMBOSS/acd
export EMBOSS_DATA
PLPLOT_LIB=$EMBOSS_PATH/share/EMBOSS
export PLPLOT_LIB

exec pepstats "$@"    
