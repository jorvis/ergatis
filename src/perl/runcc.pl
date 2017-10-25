#!/usr/bin/env perl

use strict;

print `ccomps -x $ARGV[0]`;

exit;
