package Log::LogSimpleLogger;

use strict;

my $OFF=0;
my $FATAL=1;
my $ERROR=2;
my $WARN=3;
my $INFO=4;
my $DEBUG=5;
my $ALL=6;

sub new {
    my $class = shift;

    #This singleton is an instance of Log::Simple This singleton is
    #required to have all multiple logger instances write to the same
    #log files with access to the global logger settings
    my $loggersingleton = shift;
    my $name = shift;

    my $self = bless {}, ref($class) || $class;

    $name = 'default' if(!defined $name);

    $self->{_name} = $name;

    die if(!defined $loggersingleton);
    $self->set_logger_instance($loggersingleton);
    $self->{_LOG_LEVEL} = $self->level();

    return $self;
}

sub set_logger_instance{
    my($self,$instance) = @_;
    $self->{_logsimpleobj} = $instance;
}

sub fatal{
    my($self,$msg) = @_;
    if($self->{_LOG_LEVEL} >= $FATAL || $self->{_logsimpleobj}->{_LOG_LEVEL} >= $FATAL){
	$self->{_logsimpleobj}->_output($msg,$self->{_name},$FATAL,caller());
    }
}

sub error{
    my($self,$msg) = @_;
    if($self->{_LOG_LEVEL} >= $ERROR || $self->{_logsimpleobj}->{_LOG_LEVEL} >= $ERROR){
	$self->{_logsimpleobj}->_output($msg,$self->{_name},$ERROR,caller());
    }
}

sub warn{
    my($self,$msg) = @_;
    if($self->{_LOG_LEVEL} >= $WARN || $self->{_logsimpleobj}->{_LOG_LEVEL} >= $WARN){
	$self->{_logsimpleobj}->_output($msg,$self->{_name},$WARN,caller());
    }
}

sub info{
    my($self,$msg) = @_;
    if($self->{_LOG_LEVEL} >= $INFO || $self->{_logsimpleobj}->{_LOG_LEVEL} >= $INFO){
	$self->{_logsimpleobj}->_output($msg,$self->{_name},$INFO,caller());
    }
}

sub debug{
    my($self,$msg) = @_;
    if($self->{_LOG_LEVEL} >= $DEBUG || $self->{_logsimpleobj}->{_LOG_LEVEL} >= $DEBUG){
	$self->{_logsimpleobj}->_output($msg,$self->{_name},$DEBUG,caller());
    }
}

sub logdie{
    my($self,$msg) = @_;
    $self->{_logsimpleobj}->_output($msg,$self->{_name},$FATAL,caller());
    my($package,$filename,$line,$subroutine) = caller();
    die "Died with '$msg' at $filename line $line\n";
}

sub is_fatal{
    my $self = shift;
    return ($self->{_LOG_LEVEL} >= $FATAL) || ($self->{_logsimpleobj}->{_LOG_LEVEL} >= $FATAL);
}
sub is_error{
    my $self = shift;
    return ($self->{_LOG_LEVEL} >= $ERROR) || ($self->{_logsimpleobj}->{_LOG_LEVEL} >= $ERROR);
}
sub is_warn{
    my $self = shift;
    return ($self->{_LOG_LEVEL} >= $WARN) || ($self->{_logsimpleobj}->{_LOG_LEVEL} >= $WARN);
}
sub is_info{
    my $self = shift;
    return ($self->{_LOG_LEVEL} >= $INFO) || ($self->{_logsimpleobj}->{_LOG_LEVEL} >= $INFO);
}
sub is_debug{
    my $self = shift;
    return ($self->{_LOG_LEVEL} >= $DEBUG) || ($self->{_logsimpleobj}->{_LOG_LEVEL} >= $DEBUG);
}

#
#Set log levels for named logger
#These log levels will be overridden by the global log level set in the singleton instance of
#Log::LogSimple

sub level{
    my($self,$level) = @_;
    $self->{_logsimpleobj}->level($level);
}

sub more_logging{
    my($self,$level) = @_;
    $self->{_LOG_LEVEL} += $level;
}

sub less_logging{
    my($self,$level) = @_;
    $self->{_LOG_LEVEL} -= $level;
}

1;
