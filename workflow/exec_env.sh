#!/bin/sh

# This file sets up the WF_ROOT and PATH variables, 
# which are required the command line invocation of
# workflow tools.
#
# Source this file before using any workflow tools. 

# Set WF_ROOT to the installation directory of workflow.
WF_ROOT=/usr/local/devel/ANNOTATION/workflow-2.1
export WF_ROOT

WF_TEMPLATE=/usr/local/devel/ANNOTATION/workflow-2.1/templates
export WF_TEMPLATE

WF_ROOT_INSTALL=/usr/local/devel/ANNOTATION/workflow-2.1
export WF_ROOT_INSTALL

#Set PATH
PATH=${WF_ROOT}:${WF_ROOT}/add-ons/bin:${PATH}
export PATH
