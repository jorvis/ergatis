##############################################################################
#
# Simple Makefile for workflow installation
#
# $Id$
#
# Copyright (c) 2002, The Institute for Genomic Research. All rights reserved.
#
##############################################################################

WORKFLOW_TEMPLATE_DIR	= /usr/local/devel/ANNOTATION/shared/doc/workflow
WORKFLOW_WRAPPER_DIR	= /usr/local/devel/ANNOTATION/shared/bin
WORKFLOW_RUN_DIR	= /usr/local/annotation/

all	: blastp allvsall pe cogs nucmer promer
	cat src/wfenv_bash.sh | sed 's/RUNPATHPREFIX=.*/RUNPATHPREFIX=\/$(WORKFLOW_RUN_DIR)' | sed 's/WRAPPERPATH=.*/WRAPPERPATH=\/$(WORKFLOW_WRAPPER_DIR)' | sed 's/TEMPLATEDIR=.*/TEMPLATEDIR=\/$(WORKFLOW_TEMPLATE_DIR)' > bin/wfenv_bash.sh
	cat src/set_runtime_vars.sh | sed 's/RUNPATHPREFIX=.*/RUNPATHPREFIX=\/$(WORKFLOW_RUN_DIR)' | sed 's/WRAPPERPATH=.*/WRAPPERPATH=\/$(WORKFLOW_WRAPPER_DIR)' | sed 's/TEMPLATEDIR=.*/TEMPLATEDIR=\/$(WORKFLOW_TEMPLATE_DIR)' > bin/set_runtime_vars.sh

blastp	:
	@echo "Installing blastp workflows in $(WORKFLOW_TEMPLATE_DIR)"
	cp src/blastp_template.xml $(WORKFLOW_TEMPLATE_DIR)
	cp src/blastp_multi_template.xml $(WORKFLOW_TEMPLATE_DIR)
	cp src/blastp.ini $(WORKFLOW_TEMPLATE_DIR)
	cp src/blastp_multi.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing blastp wrapper in $(WORKFLOW_WRAPPER_DIR)"
	mkdir -p bin
	cat src/run_blastp.sh | sed 's/RUNPATHPREFIX=.*/RUNPATHPREFIX=\/$(WORKFLOW_RUN_DIR)' | sed 's/WRAPPERPATH=.*/WRAPPERPATH=\/$(WORKFLOW_WRAPPER_DIR)' | sed 's/TEMPLATEDIR=.*/TEMPLATEDIR=\/$(WORKFLOW_TEMPLATE_DIR)' > bin/run_blastp.sh
	cp bin/run_blastp.sh $(WORKFLOW_WRAPPER_DIR)

allvsall	:
	@echo "Installing all_vs_all workflows in $(WORKFLOW_TEMPLATE_DIR)"
	cp src/all_vs_all_template.xml $(WORKFLOW_TEMPLATE_DIR)
	cp src/all_vs_all_multi_template.xml $(WORKFLOW_TEMPLATE_DIR)
	cp src/all_vs_all.ini $(WORKFLOW_TEMPLATE_DIR)
	cp src/all_vs_all_multi.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing all_vs_all wrapper in $(WORKFLOW_WRAPPER_DIR)"
	mkdir -p bin
	cat src/run_all_vs_all.sh | sed 's/RUNPATHPREFIX=.*/RUNPATHPREFIX=\/$(WORKFLOW_RUN_DIR)' | sed 's/WRAPPERPATH=.*/WRAPPERPATH=\/$(WORKFLOW_WRAPPER_DIR)' | sed 's/TEMPLATEDIR=.*/TEMPLATEDIR=\/$(WORKFLOW_TEMPLATE_DIR)' > bin/run_all_vs_all.sh
	cp bin/run_all_vs_all.sh $(WORKFLOW_WRAPPER_DIR)

pe	:
	@echo "Installing pe workflows in $(WORKFLOW_TEMPLATE_DIR)"
	cp src/pe_template.xml $(WORKFLOW_TEMPLATE_DIR)
	cp src/pe_multi_template.xml $(WORKFLOW_TEMPLATE_DIR)
	cp src/pe.ini $(WORKFLOW_TEMPLATE_DIR)
	cp src/pe_multi.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing pe wrapper in $(WORKFLOW_WRAPPER_DIR)"
	mkdir -p bin
	cat src/run_pe.sh | sed 's/RUNPATHPREFIX=.*/RUNPATHPREFIX=\/$(WORKFLOW_RUN_DIR)' | sed 's/WRAPPERPATH=.*/WRAPPERPATH=\/$(WORKFLOW_WRAPPER_DIR)' | sed 's/TEMPLATEDIR=.*/TEMPLATEDIR=\/$(WORKFLOW_TEMPLATE_DIR)' > bin/run_pe.sh
	cp bin/run_pe.sh $(WORKFLOW_WRAPPER_DIR)

cogs	:
	@echo "Installing cogs workflows in $(WORKFLOW_TEMPLATE_DIR)"
	cp src/cogs.xml $(WORKFLOW_TEMPLATE_DIR)
	cp src/cogs.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing cogs wrapper in $(WORKFLOW_WRAPPER_DIR)"
	mkdir -p bin
	cat src/run_cogs.sh | sed 's/RUNPATHPREFIX=.*/RUNPATHPREFIX=\/$(WORKFLOW_RUN_DIR)' | sed 's/WRAPPERPATH=.*/WRAPPERPATH=\/$(WORKFLOW_WRAPPER_DIR)' | sed 's/TEMPLATEDIR=.*/TEMPLATEDIR=\/$(WORKFLOW_TEMPLATE_DIR)' > bin/run_cogs.sh
	cp bin/run_cogs.sh $(WORKFLOW_WRAPPER_DIR)

db2bsml	:
	@echo "Installing db2bsml workflows in $(WORKFLOW_TEMPLATE_DIR)"
	cp src/db2bsml.xml $(WORKFLOW_TEMPLATE_DIR)
	cp src/db2bsml.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing db2bsml wrapper in $(WORKFLOW_WRAPPER_DIR)"
	mkdir -p bin
	cat src/run_db2bsml.sh | sed 's/RUNPATHPREFIX=.*/RUNPATHPREFIX=\/$(WORKFLOW_RUN_DIR)' | sed 's/WRAPPERPATH=.*/WRAPPERPATH=\/$(WORKFLOW_WRAPPER_DIR)' | sed 's/TEMPLATEDIR=.*/TEMPLATEDIR=\/$(WORKFLOW_TEMPLATE_DIR)' > bin/run_db2bsml.sh
	cp bin/run_db2bsml.sh $(WORKFLOW_WRAPPER_DIR)

nucmer	:
	@echo "Installing nucmer workflows in $(WORKFLOW_TEMPLATE_DIR)"
	cp src/nucmer_template.xml $(WORKFLOW_TEMPLATE_DIR)
	cp src/nucmer_multi_template.xml $(WORKFLOW_TEMPLATE_DIR)
	cp src/nucmer.ini $(WORKFLOW_TEMPLATE_DIR)
	cp src/nucmer_multi.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing nucmer wrapper in $(WORKFLOW_WRAPPER_DIR)"
	mkdir -p bin
	cat src/run_nucmer.sh | sed 's/RUNPATHPREFIX=.*/RUNPATHPREFIX=\/$(WORKFLOW_RUN_DIR)' | sed 's/WRAPPERPATH=.*/WRAPPERPATH=\/$(WORKFLOW_WRAPPER_DIR)' | sed 's/TEMPLATEDIR=.*/TEMPLATEDIR=\/$(WORKFLOW_TEMPLATE_DIR)' > bin/run_nucmer.sh
	cp bin/run_nucmer.sh $(WORKFLOW_WRAPPER_DIR)

promer	:
	@echo "Installing promer workflows in $(WORKFLOW_TEMPLATE_DIR)"
	cp src/promer_template.xml $(WORKFLOW_TEMPLATE_DIR)
	cp src/promer_multi_template.xml $(WORKFLOW_TEMPLATE_DIR)
	cp src/promer.ini $(WORKFLOW_TEMPLATE_DIR)
	cp src/promer_multi.ini $(WORKFLOW_TEMPLATE_DIR)
	@echo "Installing promer wrapper in $(WORKFLOW_WRAPPER_DIR)"
	mkdir -p bin
	cat src/run_promer.sh | sed 's/RUNPATHPREFIX=.*/RUNPATHPREFIX=\/$(WORKFLOW_RUN_DIR)' | sed 's/WRAPPERPATH=.*/WRAPPERPATH=\/$(WORKFLOW_WRAPPER_DIR)' | sed 's/TEMPLATEDIR=.*/TEMPLATEDIR=\/$(WORKFLOW_TEMPLATE_DIR)' > bin/run_promer.sh
	cp bin/run_promer.sh $(WORKFLOW_WRAPPER_DIR)
