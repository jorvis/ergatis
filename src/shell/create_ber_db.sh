#!/bin/sh

print_usage()
{
	progname=`basename $0`
	cat << END
usage: $progname -i <btab_hits> -o <nucdb_fasta_out>
	-d <db> -f <db_format> -b <bin_dir> -y <fetch_bin_dir> [-p] [-m <prot_nuc_id_map>]
END
	exit 1
}

protein=F

while getopts "i:o:d:f:m:b:p" opt
do
	case $opt in
		i) hits=$OPTARG;;
		o) out=$OPTARG;;
		d) db=$OPTARG;;
		f) db_format=$OPTARG;;
		m) id_map=$OPTARG;;
		b) bin_dir=$OPTARG;;
        y) fetch_bin_dir=$OPTARG;;
		p) protein='T';;
	esac
done

test -z $hits && echo "No btab hits provided" && print_usage
test -z $out &&	echo "No output nucleotide fasta provided" && print_usage
test -z $db && echo "No db provided" && print_usage
test -z $db_format && echo "No db format provided" && print_usage
test -z $bin_dir && echo "No bin directory provided" && print_usage
test $protein == 'F' -a -z "$id_map" && echo "No protein/nucleotide map provided" && print_usage

if [ -f $hits ] && [ ! -s $hits ]
then
	touch $out
	exit 0
fi

if [ $protein == 'F' ]
then
	prot_id=`cut -f1 $hits | head -1`
	if [ $prot_id ]
	then
		id_to_fetch=`grep "${prot_id//\./\\\.}	" $id_map | cut -f2 | sort -u`
	fi
else
	id_to_fetch=`cut -f6 $hits | sort -u | perl -pe 's/\n/ /g'`
fi

if [ "$id_to_fetch" ]
then
	$bin_dir/fetch_fasta_from_db -i "$id_to_fetch" -d $db -p $protein -f $db_format -o $out -b $fetch_bin_dir
	ec=$?
	if [ $ec -ne 0 ]
	then
		echo "Error fetching sequences"
		exit $ec
	fi
else
	exit 1
fi
