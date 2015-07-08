package Ergatis::IdGenerator::Config;

## This config file defines the module that will be used for id generation by ergatis
## The specified class should at least implement all methods defined in Ergatis::IdGenerator

use Ergatis::IdGenerator::IGSIdGenerator;
our $class = 'Ergatis::IdGenerator::IGSIdGenerator';

1;
