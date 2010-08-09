package MG::Common;


use vars qw(@ISA @EXPORT $VERSION);
use Exporter;

@ISA = qw(Exporter);
@EXPORT  = qw(&initialize);

sub initialize {
#directories
    $global{local_bin} = '/usr/local/bin/';
    $global{tmp} = '/tmp';
    $global{scratch} = '/local/scratch/';
    $global{metagene_bin} = '/export/seqprg/bin/metagene';

#blast parameters
    $global{evalue} = '1e-5';
    $global{num_align} = 150;
    $global{sffbin} = '/usr/local/packages/offInstrumentApps_2.0.01.12-64/bin/';
    return %global;
}
