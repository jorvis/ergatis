##############################################################################
#
# Simple Makefile for workflow installation
#
# $Id$
#
# Copyright (c) 2002, The Institute for Genomic Research. All rights reserved.
#
##############################################################################

#WORKFLOW_TEMPLATE_DIR	= /usr/local/devel/ANNOTATION/shared/doc/workflow
#WORKFLOW_WRAPPER_DIR	= /usr/local/devel/ANNOTATION/shared/bin
#WORKFLOW_RUN_DIR	= /usr/local/annotation/

WORKFLOW_TEMPLATE_DIR	= /tmp/workflow/workflow
WORKFLOW_WRAPPER_DIR	= /tmp/workflow/bin
WORKFLOW_RUN_DIR	= /tmp/workflow/run

all	: install

dirs	:
	@mkdir -p $(WORKFLOW_TEMPLATE_DIR)
	@mkdir -p $(WORKFLOW_WRAPPER_DIR)
	@mkdir -p $(WORKFLOW_RUN_DIR)

install	: dirs db2bsml blastp allvsall pe cogs nucmer promer
	@cat src/wfenv_bash.sh | sed 's?RUNPATHPREFIX=.*?RUNPATHPREFIX=$(WORKFLOW_RUN_DIR)?' | sed 's?WRAPPERPATH=.*?WRAPPERPATH=$(WORKFLOW_WRAPPER_DIR)?' | sed 's?TEMPLATEDIR=.*?TEMPLATEDIR=$(WORKFLOW_TEMPLATE_DIR)?' > bin/wfenv_bash.sh
	@cat src/set_runtime_vars.sh | sed 's?RUNPATHPREFIX=.*?RUNPATHPREFIX=$(WORKFLOW_RUN_DIR)?' | sed 's?WRAPPERPATH=.*?WRAPPERPATH=$(WORKFLOW_WRAPPER_DIR)?' | sed 's?TEMPLATEDIR=.*?TEMPLATEDIR=$(WORKFLOW_TEMPLATE_DIR)?' > bin/set_runtime_vars.sh
	@install bin/wfenv_bash.sh $(WORKFLOW_WRAPPER_DIR)
	@install bin/set_runtime_vars.sh $(WORKFLOW_WRAPPER_DIR)

blastp	:
	@echo "Installing blastp workflows in $(WORKFLOW_TEMPLATE_DIR)"
	@install src/blastp_template.xml $(WORKFLOW_TEMPLATE_DIR)
	@install src/blastp_multi_template.xml $(WORKFLOW_TEMPLATE_DIR)
	@install src/blastp.ini $(WORKFLOW_TEMPLATE_DIR)
	@install src/blastp_multi.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing blastp wrapper in $(WORKFLOW_WRAPPER_DIR)"
	@mkdir -p bin
	@cat src/run_blastp.sh | sed 's?RUNPATHPREFIX=.*?RUNPATHPREFIX=$(WORKFLOW_RUN_DIR)?' | sed 's?WRAPPERPATH=.*?WRAPPERPATH=$(WORKFLOW_WRAPPER_DIR)?' | sed 's?TEMPLATEDIR=.*?TEMPLATEDIR=$(WORKFLOW_TEMPLATE_DIR)?' > bin/run_blastp.sh
	@install bin/run_blastp.sh $(WORKFLOW_WRAPPER_DIR)

allvsall	:
	@echo "Installing all_vs_all workflows in $(WORKFLOW_TEMPLATE_DIR)"
	@install src/all_vs_all_template.xml $(WORKFLOW_TEMPLATE_DIR)
	@install src/all_vs_all_multi_template.xml $(WORKFLOW_TEMPLATE_DIR)
	@install src/all_vs_all.ini $(WORKFLOW_TEMPLATE_DIR)
	@install src/all_vs_all_multi.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing all_vs_all wrapper in $(WORKFLOW_WRAPPER_DIR)"
	@mkdir -p bin
	@cat src/run_all_vs_all.sh | sed 's?RUNPATHPREFIX=.*?RUNPATHPREFIX=$(WORKFLOW_RUN_DIR)?' | sed 's?WRAPPERPATH=.*?WRAPPERPATH=$(WORKFLOW_WRAPPER_DIR)?' | sed 's?TEMPLATEDIR=.*?TEMPLATEDIR=$(WORKFLOW_TEMPLATE_DIR)?' > bin/run_all_vs_all.sh
	@install bin/run_all_vs_all.sh $(WORKFLOW_WRAPPER_DIR)

pe	:
	@echo "Installing pe workflows in $(WORKFLOW_TEMPLATE_DIR)"
	@install src/pe_template.xml $(WORKFLOW_TEMPLATE_DIR)
	@install src/pe_multi_template.xml $(WORKFLOW_TEMPLATE_DIR)
	@install src/pe.ini $(WORKFLOW_TEMPLATE_DIR)
	@install src/pe_multi.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing pe wrapper in $(WORKFLOW_WRAPPER_DIR)"
	@mkdir -p bin
	@cat src/run_pe.sh | sed 's?RUNPATHPREFIX=.*?RUNPATHPREFIX=$(WORKFLOW_RUN_DIR)?' | sed 's?WRAPPERPATH=.*?WRAPPERPATH=$(WORKFLOW_WRAPPER_DIR)?' | sed 's?TEMPLATEDIR=.*?TEMPLATEDIR=$(WORKFLOW_TEMPLATE_DIR)?' > bin/run_pe.sh
	@install bin/run_pe.sh $(WORKFLOW_WRAPPER_DIR)

cogs	:
	@echo "Installing cogs workflows in $(WORKFLOW_TEMPLATE_DIR)"
	@install src/cogs_template.xml $(WORKFLOW_TEMPLATE_DIR)
	@install src/cogs.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing cogs wrapper in $(WORKFLOW_WRAPPER_DIR)"
	@mkdir -p bin
	@cat src/run_cogs.sh | sed 's?RUNPATHPREFIX=.*?RUNPATHPREFIX=$(WORKFLOW_RUN_DIR)?' | sed 's?WRAPPERPATH=.*?WRAPPERPATH=$(WORKFLOW_WRAPPER_DIR)?' | sed 's?TEMPLATEDIR=.*?TEMPLATEDIR=$(WORKFLOW_TEMPLATE_DIR)?' > bin/run_cogs.sh
	@install bin/run_cogs.sh $(WORKFLOW_WRAPPER_DIR)

db2bsml	:
	@echo "Installing db2bsml workflows in $(WORKFLOW_TEMPLATE_DIR)"
	@install src/db2bsml_template.xml $(WORKFLOW_TEMPLATE_DIR)
	@install src/db2bsml.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing db2bsml wrapper in $(WORKFLOW_WRAPPER_DIR)"
	@mkdir -p bin
	@cat src/run_db2bsml.sh | sed 's?RUNPATHPREFIX=.*?RUNPATHPREFIX=$(WORKFLOW_RUN_DIR)?' | sed 's?WRAPPERPATH=.*?WRAPPERPATH=$(WORKFLOW_WRAPPER_DIR)?' | sed 's?TEMPLATEDIR=.*?TEMPLATEDIR=$(WORKFLOW_TEMPLATE_DIR)?' > bin/run_db2bsml.sh
	@install bin/run_db2bsml.sh $(WORKFLOW_WRAPPER_DIR)

nucmer	:
	@echo "Installing nucmer workflows in $(WORKFLOW_TEMPLATE_DIR)"
	@install src/nucmer_template.xml $(WORKFLOW_TEMPLATE_DIR)
	@install src/nucmer_multi_template.xml $(WORKFLOW_TEMPLATE_DIR)
	@install src/nucmer.ini $(WORKFLOW_TEMPLATE_DIR)
	@install src/nucmer_multi.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing nucmer wrapper in $(WORKFLOW_WRAPPER_DIR)"
	@mkdir -p bin
	@cat src/run_nucmer.sh | sed 's?RUNPATHPREFIX=.*?RUNPATHPREFIX=$(WORKFLOW_RUN_DIR)?' | sed 's?WRAPPERPATH=.*?WRAPPERPATH=$(WORKFLOW_WRAPPER_DIR)?' | sed 's?TEMPLATEDIR=.*?TEMPLATEDIR=$(WORKFLOW_TEMPLATE_DIR)?' > bin/run_nucmer.sh
	@install bin/run_nucmer.sh $(WORKFLOW_WRAPPER_DIR)

promer	:
	@echo "Installing promer workflows in $(WORKFLOW_TEMPLATE_DIR)"
	@install src/promer_template.xml $(WORKFLOW_TEMPLATE_DIR)
	@install src/promer_multi_template.xml $(WORKFLOW_TEMPLATE_DIR)
	@install src/promer.ini $(WORKFLOW_TEMPLATE_DIR)
	@install src/promer_multi.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing promer wrapper in $(WORKFLOW_WRAPPER_DIR)"
	@mkdir -p bin
	@cat src/run_promer.sh | sed 's?RUNPATHPREFIX=.*?RUNPATHPREFIX=$(WORKFLOW_RUN_DIR)?' | sed 's?WRAPPERPATH=.*?WRAPPERPATH=$(WORKFLOW_WRAPPER_DIR)?' | sed 's?TEMPLATEDIR=.*?TEMPLATEDIR=$(WORKFLOW_TEMPLATE_DIR)?' > bin/run_promer.sh
	@install bin/run_promer.sh $(WORKFLOW_WRAPPER_DIR)
