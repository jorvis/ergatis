package PFunc::TIGRRolesAnnotation;

use strict;
use warnings;
use File::OpenFile qw( open_file );
use Data::Dumper;

#CONSTANTS
my $uncategorized_tigr_role = 185;

sub new {
    my ($class, %args) = @_;
    my $self = {};
    bless( $self, $class );
    
    if( $args{'tigr_roles_db_dir'} ) {
        $self->{'_keyword_file'} = $args{'tigr_roles_db_dir'}."/tigr_roles_keywords.txt";
        die("Could not find keywords file: $self->{'_keyword_file'}") 
            unless( -e $self->{'_keyword_file'} );
    } else {
        die("Argument tigr_roles_db_dir is required for creation of a TIGRRolesAnnotation object");
    }
    return $self;
}

sub assign_tigr_roles_by_keyword {
    my ($self, $common_name, $cur_roles, $flag) = @_;
    my $search_keywords = 0;
    my $new_tigr_roles;

    print "Common name: $common_name\n" if( $flag );

    if( !defined( $cur_roles ) || @{$cur_roles} == 0 || $cur_roles->[0] == $uncategorized_tigr_role ) {
        $search_keywords = 1;
    }

    if( $search_keywords ) {
        $self->{'_tigr_role_lookup'} = &_parse_keyword_file( $self->{'_keyword_file'} )
            unless( exists( $self->{'_tigr_role_lookup'} ) );
        
        foreach my $map ( reverse( @{$self->{'_tigr_role_lookup'}} ) ) {
            my $keyword = $map->{'keyword'};

            if( !defined( $keyword ) ) {
                die("keyword is not defined");
            }

            if( !defined( $common_name ) ) {
                die("Don't have a common name");
            }

            if( $common_name =~ /$keyword/i ) {
                $new_tigr_roles = [$map->{'role_id'}];
            }
            
        }

        #if at this point we still don't have any TIGR roles
        if( !defined($new_tigr_roles) || @{$new_tigr_roles} == 0 ) {
            $new_tigr_roles = [$uncategorized_tigr_role];
        }
    }

    print "returning ".join(" ", @{$new_tigr_roles})."\n" if( $flag && defined( $new_tigr_roles ) );
    return $new_tigr_roles;
    
}


sub _parse_keyword_file {
    my ($keywords_file) = @_;
    my $retval = [];
    
    my $index = 0;
    my $in = &open_file( $keywords_file, 'in' );
    
    while( my $line = <$in> ) {
        chomp( $line );
        my @cols = split(/\s+/, $line );
        my $role_id = pop @cols;
        my $keyword = join( " ", @cols );
        die("keyword not defined [$line]") unless( $keyword );

        push(@{$retval}, {
            'role_id' => $role_id,
            'keyword' => $keyword
            });
    }
    
    return $retval;
}
