#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Basename;
use JSON;
use POSIX qw(strftime);
$| = 1;

my %options = ();
my $results = GetOptions (\%options,
              'input_list|i=s',
              'tags_to_upload|ttu=s',
              'output_directory|o=s',
              'data_output_directory|do=s',
              'credential|c=s',
              'user|u=s',
              'key|k=s',
              'host|s=s',
              'decrypt_script|d=s',
              'skip_transfer|st=s',
              'child_config|cc=s',
              'child_config_params|ccp=s',
              'pipeline_parent|pp=s',
              'num_retries|nr=s',
              'input_param_name|ipn=s',
              'wait_time|wt=s',
              'tags_to_download|ttd=s',
              'max_pipelines|mp=s',
              'help|h');

my $pipe_tasks = [];
my $input_lists = undef;
my $finished = 0;
my @tags_to_upload = split(/\,/, $options{tags_to_upload});
$options{output_directory} =~ s/\/$//;
my $data_directory = "$options{output_directory}/data";
my $task_file = $options{output_directory}."/pipeline_tasks.txt";
my $wait = $options{wait_time} ? $options{wait_time} : 100;
my $max_pipes = $options{max_pipelines} ? $options{max_pipelines} : 50;
my $max_pipes_file = "$options{output_directory}/throttle.txt";
my $resize_interval = 3600;
my $last_resize = time;
my $num_retries = $options{num_retries} ? $options{num_retries} : 2;
    
my $num_running = 0;

while(!$finished) {

    # Read in the pipeline tasks if we haven't already
    if(-e $task_file && !(@$pipe_tasks)) {
        print &get_time()." - LOG: Loading task file\n";
        $pipe_tasks = &load_state($task_file);
        map{
            if(!$_->{completed}) {
                $num_running++;
            }
        }@$pipe_tasks;
    }
    
    # Load the input list of lists
    if(not defined $input_lists) {
        print &get_time()." - LOG: Loading Input Lists\n";
        &load_input_lists();
    }
    if(-e $max_pipes_file) {
        print &get_time()." - LOG: Found throttle file\n";
        my $num = `cat $max_pipes_file`;
        $num =~ s/[\D]//g;
        $max_pipes = $num;
        print &get_time()." - LOG: now throttling to $max_pipes pipelines\n";
    }
    # Check if any tasks have finished and download the data if they have
    $finished = &check_tasks($pipe_tasks);
    
    # Submit the next task
    if(@$input_lists && ($num_running < $max_pipes)) {
        &submit_pipeline(shift(@$input_lists));
    }
    
    # Check that all clusters are the right size (avoids losing all execs)
    if(time - $last_resize >= $resize_interval) {
        # Check if any tasks have finished and download the data if they have
        # Need to do this here as well. Don't want to resize if we are finished.
        $finished = &check_tasks($pipe_tasks);
        
        &resize_clusters();
        $last_resize = time;
    }
    print &get_time()." - LOG: Dumping Tasks at the end of an iteration\n";    
    &dump_tasks();
    print &get_time()." - LOG: Finished with an iteration, waiting $wait seconds\n";
    sleep $wait;
}

print &get_time()." - LOG: Finished!\n";

exit(0);

sub get_time {
    return strftime "%a %b %e %H:%M:%S %Y", gmtime;
}

sub submit_pipeline {
    my $input = shift;
    my $task = shift;
    my $fname = fileparse($input,qr/\.[^.]*$/);
    print &get_time()." - LOG: Going to submit $fname since $num_running < $max_pipes\n";
    if(!defined($task) && $num_running >= $max_pipes) {
        return;
    }
    
    # Create an output directory for the temp files for this subgroup
    my $output_dir = "$options{output_directory}/$fname";
    if(! -d $output_dir) {
        my $cmd = "mkdir $output_dir";
       `$cmd`;
        if($?) {
            die &get_time()." - FAILURE: Failed running:\n\t$cmd\n$?";
        }
    }
    else {
        print &get_time()." - WARNING: output directory $output_dir already exists\n";
    }
    
    # Start new cluster
    my $cmd = "vp-start-cluster --cluster=$fname\_cluster --num-exec=0 --cred=$options{credential} -t";
   `$cmd`;
    if($?) {
        &print_error($cmd,$input);
        return;
    }
    
    # Copy files and decrypt files (if necessary)
    my $localfl = "$output_dir/$fname.list";
    print &get_time()." - LOG: Local file list will be $localfl\n";
    &copy_files($input,$localfl);
    
    # Create tag
    print &get_time()." - LOG: Tagging $fname\n";
    &create_tag($fname,$localfl);
    
    # Check cluster status, wait till it's up to continue
    print &get_time()." - LOG: Waiting for cluster $fname\_cluster to come up\n";
    my $cmd = "vp-start-cluster --cluster=$fname\_cluster --num-exec=0 --cred=$options{credential}";
   `$cmd`;
    if($?) {     
        &print_error($cmd,$input);
        # WORST THING EVER, submitting a bogus cluster to take the place of a good cluster.
        print  &get_time()." - LOG: Submitting bogus cluster $fname\_bogus to try to recover\n";
        my $cmd = "vp-start-cluster --cluster=$fname\_bogus --num-exec=0 --cred=$options{credential} -t";
        `$cmd`;
        if($?) {
            &print_error($cmd);
        }
        return;
    }
    print &get_time()." - LOG: Cluster $fname\_cluster came up\n";
    
    # Configure queues on exec nodes
    print &get_time()." - LOG: Getting the master hostname for cluster $fname\_cluster\n";
    my $cmd = "vp-describe-cluster --cluster=$fname\_cluster";
    my @ret = `$cmd`;
    if($?) {
        &print_error($cmd,$input);
        return;
    }
    my $master = '';
    map {if(/MASTER\s+\S+\s+(\S+)\s+/) {$master = $1;}}@ret;
    print &get_time()." - LOG: Setting the exec queues to 1 slot and turning off autoshutdown for master host $master on $fname\_cluster\n";
    my $cmd = "/opt/ergatis/global_saved_templates/clovr_lgt_wrapper/set_queue.sh $master";
    `$cmd`;
    if($?) {
        &print_error($cmd,$input);
        return;
    }    
    
    # Add instances
    my $num_lines = `wc -l $localfl`;
    $num_lines =~ s/(\d+)\s+.*/$1/;
    chomp $num_lines;
    my $num_execs = $num_lines;
    print &get_time()." - LOG: Resizing cluster $fname\_cluster with $num_execs execs\n";
    my $cmd = "vp-run-metrics -c cluster.CLUSTER_NAME=$fname\_cluster -c pipeline.EXEC_WANT_INSTANCES=$num_execs resize-cluster";
   `$cmd`;
    if($?) {
        &print_error($cmd,$input);
        return;
    }    
    # Upload tags
    print &get_time()." - LOG: Uploading several tags to $fname\_cluster\n";
    &upload_tags([$fname,@tags_to_upload],"$fname\_cluster");
    # Write configuration
    
    print &get_time()." - LOG: Writing config options to $output_dir/$fname\_pipeline.config\n";
    my $ipn = $options{input_param_name} ? $options{input_param_name} : 'input.INPUT_TAG';
    my $cmd = "vp-describe-protocols --config-from-protocol=$options{child_config} -c $ipn=$fname -c pipeline.PIPELINE_WRAPPER_NAME=$fname $options{child_config_params} > $output_dir/$fname\_pipeline.config";
    print &get_time()." - LOG: running $cmd\n";
    `$cmd`;
    if($?) {
        &print_error($cmd,$input);
        return;
    }    
    # Run Pipeline
    print &get_time()." - LOG: Running pipeline with $output_dir/$fname\_pipeline.config on $fname\_cluster\n";

    my $cmd = "vp-run-pipeline --pipeline-config=$output_dir/$fname\_pipeline.config --pipeline-queue=pipeline.q --bare --cluster=$fname\_cluster -t --overwrite";
    if($options{pipeline_parent}) {
        $cmd .= " --pipeline-parent=$options{pipeline_parent}";
    }
    my $tid = `$cmd`;
    chomp $tid;
    if($?) {
        &print_error($cmd,$input);
        return;
    }      

    my $downloads = [];
    map{
        push(@$downloads,"$fname\_$_");
    } split(/,/,$options{tags_to_download});

    my $newtask = {
            input_list => $input,
            task_id => $tid,
            completed => 0,
            num_retries => $num_retries,
            num_task_retries => $num_retries,
            cluster => "$fname\_cluster",
            download_tags => $downloads,
            num_execs => $num_execs
        };
    # If a task was passed in then we are resubmitting    
    if(!$task) {
        # Create pipeline task entry and write out pipeline tasks
        push(@$pipe_tasks,$newtask);
        $num_running++;
    }
    else {
        print &get_time()." - LOG: Ran $fname on $fname\_cluster with task_id: $tid Old task_id: $task->{task_id}\n";
        $newtask->{num_retries} = $task->{num_retries};
        $newtask->{num_task_retries} = $task->{num_retries};        
        %$task = %$newtask;
    }
    &dump_tasks();    
}

sub print_error {
    my $cmd = shift;
    my $file = shift;
    print &get_time()." - FAILURE: Failed running:\n$cmd\n$?";    
    if($file) {
        print &get_time()." - LOG: Putting $file back on the list\n";
        unshift(@$input_lists, $file);        
    }
}

sub load_input_lists {
    my $submitted_lists = {};
    map {$submitted_lists->{$_->{input_list}} = 1}@$pipe_tasks;
    $input_lists = [];
    open IN, "<$options{input_list}" or die "Unable to open $options{input_list}\n";
    while(<IN>) {
        chomp;
        if(!$submitted_lists->{$_}) {
            print &get_time()." - LOG: Need to process $_\n";
            push(@$input_lists, $_);
        }
    }
}

sub resize_clusters {

    print &get_time()." - LOG: Checking Cluster sizes\n";
    foreach my $task (@{$pipe_tasks}) {
        if(!$task->{completed}) {
            my $cmd = "vp-run-metrics -c cluster.CLUSTER_NAME=$task->{cluster} -c pipeline.EXEC_WANT_INSTANCES=$task->{num_execs} resize-cluster";
           `$cmd`;
            if($?) {
               print &get_time()." - FAILURE: Failed running:\n$cmd\n$?";
            }
        }
    }
}

sub load_state {
    my $file = shift;
    
    open IN, "<$file" or die "Unable to open $file\n";
    my @text = <IN>;
    close IN;
    my $str = join('',@text);
    return decode_json($str);
}

sub check_tasks {
    my $tasks = shift;
    my $num_complete = 0;
    my $finished = 0;
    foreach my $task (@$tasks) {
    
        # If we have completed a task, harvest the data
        # and update the task list. If the task failed then
        # mark it complete so we can finish
        if(!$task->{completed}) {
            my $complete = &check_task_complete($task);
            
            if($complete) {
                if($complete ne 'failed') {
                    print "$task->{task_id} complete on cluster $task->{cluster}\n";
                    print "Harvesting data from $task->{cluster}\n";
                    &harvest_data($task);
                }
                &terminate_cluster($task->{cluster});
                $task->{completed} = 1;
                &dump_tasks();
                $num_running--;
            }
        }
        if($task->{completed}) {
            $num_complete++;
        }
    }
    if($num_complete == @$tasks && scalar @$input_lists ==0) {
        print &get_time()." - LOG: It looks like we're done, $num_complete complete of ".(scalar @$tasks)." and ".(scalar @$input_lists)." inputs left\n"; 
        $finished = 1;
    }
    return $finished;
}

sub dump_tasks {
    open OUT, ">$task_file" or die "Unable to open $task_file\n";
    print OUT to_json($pipe_tasks,{pretty => 1});
    close OUT;
}

sub terminate_cluster {
    my $cluster = shift;
    
    print &get_time()." - LOG: Checking status of tasks on $cluster\n";
    my $cmd = "vp-describe-task --name=$cluster --block";
    `$cmd`;
    my $ret = `$cmd`;
    if($?) {
        print &get_time()." - WARNING: Task status not available for $cluster, will terminate anyway\n$?";
    }
    print &get_time()." - LOG: Terminating cluster $cluster\n";
    my $cmd = "vp-terminate-cluster --cluster=$cluster";
    my $ret = `$cmd`;
    if($?) {
        print &get_time()." - WARNING: Failed to terminate cluster $cluster, will move on anyway\n$?";
    }
}

sub clear_data {

    if(-d $data_directory) {
        print &get_time()." - LOG: Clearing out the data directory\n";
        my $cmd = "rm $data_directory/*";
        `$cmd`;
        if($?) {
            print &get_time()." - WARNING: Had a problem clearing the output directory\n$?\n";
        }
    }
    if(-d "$data_directory/decrypt") {
        print &get_time()." - LOG: Clearing the decrypt directory\n";
        my $cmd = "rm $data_directory/decrypt/*";
        `$cmd`;
        if($?) {
            print &get_time()." - WARNING: Had a problem clearing the decrypt directory\n$?\n";
        }
    }
}

sub check_task_complete {
    my $task = shift;
    my $tid = $task->{task_id};
    my $cluster = $task->{cluster};
    my $complete = 0;
    print &get_time()." - LOG: Checking status of task $tid running on $cluster\n";
    my $cmd = "vp-describe-task $tid --name=$cluster";
    my $ret = `$cmd`;
    if($?) {
        print &get_time()." - FAILURE: Failed to get info about task $tid on $cluster:\n$cmd\n$?";
        if($task->{num_task_retries} > 0) {
            print &get_time()." - LOG: Couldn't get info for $tid on $cluster will try again at the next iteration. $task->{num_task_retries} retries remain.\n";
            $task->{num_task_retries}--;
        }
        elsif($task->{num_retries} > 0) {
            $task->{num_retries}--;
            &terminate_cluster($task->{cluster});
            print &get_time()." - LOG: Retrying $tid on $cluster $task->{num_retries} remain\n";
            &submit_pipeline($task->{input_list},$task);
        }
        else {
            $complete = 'failed';
        }
    }
    else {
        $task->{num_task_retries} = $num_retries;
    }
    
    if($ret =~ /State: complete/) {
        $complete = 1;
    }
    elsif($ret =~ /State: failed/) {
        if($task->{num_retries} > 0) {
            $task->{num_retries}--;
            print &get_time()." - LOG: Retrying $tid on $cluster $task->{num_retries} remain\n";
            &submit_pipeline($task->{input_list},$task);
        }
        else {
            $complete = 'failed';
        }
    }
    return $complete;
}

sub harvest_data {
    my $task = shift;
    print &get_time()." - LOG: Harvesting data from $task->{cluster}\n";
    &clear_data();
    foreach my $tag (@{$task->{download_tags}}) {
        # First step is transfering the dataset back to the local node.
        print &get_time()." - LOG: transfering $tag from $task->{cluster}\n";
        my $cmd = "vp-transfer-dataset --src-cluster=$task->{cluster} --dst-cluster=local --cluster=local --block --tag-name=$tag";
        `$cmd`;
        if($?) {
            $task->{num_retries}--;
            print &get_time()." - FAILURE: Failed running:\n$cmd\n$?";
            next;
        }
        my @files = `vp-describe-dataset --tag-name=$tag|grep FILE|cut -f2`;
        foreach my $file (@files) {
            chomp $file;
            print &get_time()." - LOG: Transfering file $file to $options{user}\@$options{host}:$options{data_output_directory}/$task->{cluster}/\n";
            &transfer_file({
                src => $file,
                dst => "$options{user}\@$options{host}:$options{data_output_directory}/$tag/"
            });
            print &get_time()." - LOG: Removing $file from the local master.\n";
            my $cmd = "rm $file";
            `$cmd`;
            if($?) {
               print &get_time()." - FAILURE: Couldn't remove $file\n$cmd\n$?";
            }            
        }
    }
}

sub transfer_file {
    my $config = shift;

    my $cmd = "rsync -rlptD -e \"ssh -o PasswordAuthentication=no -o ConnectTimeout=30 -o StrictHostKeyChecking=no -o ServerAliveInterval=30 -o UserKnownHostsFile=/dev/null -i $options{key}\" $config->{src} $config->{dst}";

    print &get_time()." - LOG: Running:\n$cmd\n";
    `$cmd`;

    if($?) {
       print &get_time()." - FAILURE: Failed running:\n$cmd\n$?";
    }

}

sub copy_files {
    my $file_list = shift;
    my $output_list = shift;
    if(! -d "$data_directory") {
        print &get_time()." - LOG: Making temporary data directory in $data_directory\n";
        `mkdir $data_directory`;
    }
    if(! -d "$data_directory/decrypt") {
        print &get_time()." - LOG: Making temporary decrypt directory in $data_directory\n";
        `mkdir $data_directory/decrypt`;
    }
    print &get_time()." - LOG: Copying a new list of files $file_list from the data server\n";
    &clear_data();
    open IN, "<$file_list" or die "unable to open list\n";
    open OUT, ">$output_list" or die "unable to open output list\n";
    while(<IN>) {
        chomp;
        if(!$options{skip_transfer}) {
            print "Copying $_\n";
            &transfer_file({
                src => "$options{user}\@$options{host}:$_",
                dst => "$data_directory",
                key => $options{key}
            });
        }
        
        my $fname = basename($_);
        my $decrypt_fname = $fname;

        if($options{decrypt_script}) {
            print "Decrypting $data_directory/$fname\n";
            my $cmd = "$options{decrypt_script} $data_directory/$fname -out-dir $data_directory/decrypt -remove-encrypted";
            print "Running:\n $cmd\n";
            my $res = `$cmd`; 
            print $res;
            if($res !~ /File validation is skipped/ && $?) {
                die "$?\n";
            }
            my $dfname = `ls $data_directory/decrypt/`;
        	chomp $dfname;
            my $cmd = "mv $data_directory/decrypt/$dfname $data_directory";
            print "Running $cmd\n";
            my $res = `$cmd`;
            if( $?) {
                die "$?\n";
            }
            $decrypt_fname = basename $dfname;
        }
        elsif(! -e "$data_directory/$decrypt_fname") {
            my $fname = basename($decrypt_fname,'.ncbi_enc');
            $decrypt_fname = $fname;
        }
        $options{output_dir} =~ s/\/$//;
        print OUT "$data_directory/$decrypt_fname\n";
    }
    close OUT;
}

sub create_tag {
    my $tag_name = shift;
    my $filelist = shift;
    
    print &get_time()." - LOG: Tagging $filelist as $tag_name\n";
    my $cmd = "cat $filelist | vp-add-dataset --tag-name=$tag_name -o --stdin";
    `$cmd`;
    if( $?) {
       die &get_time()." - FAILURE: Failed running:\n$cmd\n$?";
    }
}

sub upload_tags {
    my $tags = shift;
    my $dest = shift;
    foreach my $tag (@$tags) {
        print &get_time()." - LOG: Transfering $tag to $dest\n";
        map {
#            my $cmd = "/opt/ergatis/global_saved_templates/clovr_lgt_wrapper/upload_list.pl $tag $dest";
            my $cmd = "vp-transfer-dataset --tag-name=$_ --dst-cluster=$dest --block --expand";
            `$cmd`;
            if( $?) {
                die &get_time()." - FAILURE: Failed running:\n$cmd\n$?";
            }
        }split(/,/, $tag)
    }
}
