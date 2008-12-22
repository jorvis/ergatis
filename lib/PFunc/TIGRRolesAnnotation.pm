package PFunc::TIGRRolesAnnotation;

use strict;
use warnings;
use File::OpenFile qw( open_file );
use Exporter 'import';
use Data::Dumper;
use vars qw(@EXPORT_OK);
@EXPORT_OK = qw( assign_tigr_roles_by_keyword );

my $uncategorized_tigr_role = 185;
my $keyword_file = "/home/kgalens/data/tigr_roles_keywords.txt";

sub assign_tigr_roles_by_keyword {
    my ($common_name, $cur_roles) = @_;
    my $search_keywords = 0;
    my $new_tigr_roles;

    if( !defined( $cur_roles ) || @{$cur_roles} == 0 || $cur_roles->[0] == 185 ) {
        $search_keywords = 1;
    }

    if( $search_keywords ) {
        my $keywords_tigr_role_lookup = &_parse_keyword_file( $keyword_file );
        
        foreach my $keyword ( keys %{$keywords_tigr_role_lookup} ) {

            if( !defined( $keyword ) ) {
                die("keyword is not defined");
            }

            if( !defined( $common_name ) ) {
                die("Don't have a common name");
            }

            if( $common_name =~ /$keyword/ ) {
                $new_tigr_roles = [$keywords_tigr_role_lookup->{ $keyword }];
            }
            
        }
    }

    return $new_tigr_roles;
    
}

sub _parse_keyword_file {
    my ($keywords_file) = @_;
    my $retval = {};
    
    my $in = &open_file( $keywords_file, 'in' );
    
    while( my $line = <$in> ) {
        chomp( $line );
        my @cols = split(/\s+/, $line );
        my $role_id = pop @cols;
        my $keyword = join( " ", @cols );
        die("keyword not defined [$line]") unless( $keyword );

        $retval->{$keyword} = $role_id;
    }
    
    return $retval;
}
