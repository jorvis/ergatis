#!/bin/sh
donefile=$1
shift
/usr/local/devel/ANNOTATION/shared/bin/linux/peffect $*
echo 'done' > $donefile
