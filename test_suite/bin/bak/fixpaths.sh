#!/bin/sh
find $1 -type f -exec chmod +w {} \;
REPPATH=$1
export REPPATH
find $1 -type f -exec perl -pi -e 's/\/usr\/local\/devel\/ANNOTATION\/ard\/[\w-]+\/lib/$ENV{"REPPATH"}\/lib/g' {} \;
