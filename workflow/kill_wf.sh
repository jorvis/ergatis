#!/bin/sh

ps | grep java | perl -ne '($id) = ($_ =~ /\s*(\d+)/);`kill -9 $id`;'
