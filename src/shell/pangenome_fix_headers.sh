#!/bin/sh
if [ -z "$SCHEMA_DOCS_DIR" ]
then
    SCHEMA_DOCS_DIR=
export SCHEMA_DOCS_DIR
fi
if [ -z "$WORKFLOW_WRAPPERS_DIR" ]
then
    WORKFLOW_WRAPPERS_DIR=/usr/local/projects/ergatis/package-nightly//bin
export WORKFLOW_WRAPPERS_DIR
fi
if [ -z "$WORKFLOW_DOCS_DIR" ]
then
    WORKFLOW_DOCS_DIR=
export WORKFLOW_DOCS_DIR
fi

umask 0000

unset LD_LIBRARY_PATH

LANG=C
export LANG
LC_ALL=C
export LC_ALL

while [[ $# -ge 1 ]]
do
i="$1"
arg=$(echo $i | cut -f1 -d "=")
val=$(echo $i | cut -f2 -d "=")

case $arg in
    --input_file)
    input_file="$val"
    ;;
    --output_file)
    output_file="$val"
    ;;
    --mapping_file)
    mapping_file="$val"
    ;;
    --python_exec)
    python_exec="$val"
    ;;
	--script_bin)
	script_bin="$val"
	;;
esac
shift
done

# This script needs to run in Python 2.7 or later but I'd like to keep a default python path in case
if [ -z "$python_exec" ]; then
	python_exec="/usr/bin/python"
fi

cmd="$python_exec ${script_bin}/pangenome_fix_headers.py --input_file=${input_file} --output_file=${output_file} --mapping_file=${mapping_file}"

echo "$cmd"
$cmd

exit 0
