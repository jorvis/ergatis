#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser set_message);

use Tree::DAG_Node;

#edit fields: name,state,type,retryCount,retryAttempts,timeOut,param.value

