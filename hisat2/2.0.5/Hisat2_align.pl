#!/usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my (@file_query, $database_path, $user_database_path, $annotation_path, 
$user_annotation_path, $file_names, $root_names, @file_query2, $file_type, $n_cores);


GetOptions( "file_query=s"      => \@file_query,
	          "file_query2=s"     => \@file_query2,
	          "database=s"        => \$database_path,
	          "user_database=s"   => \$user_database_path,
            "annotation=s"      => \$annotation_path,
            "user_annotation=s" => \$user_annotation_path,
	          "file_names=s"      => \$file_names,
	          "root_names=s"      => \$root_names,
	          "file_type=s"       => \$file_type,
            "n_cores=s"         => \$n_cores,
          );

# sanity check for input data
if (@file_query2) {
    @file_query && @file_query2 || die "Error: At least one file for each paired-end is required\n"; 
    @file_query == @file_query2 || die "Error: Unequal number of files for paired ends\n";
}

if (!($user_database_path || $database_path)) {
    die "No reference genome was supplied\n";
}
if (@file_query < 1) {
    die "No FASTQ files were supplied\n";
}

# Sanity check for input ref. genome
unless ($database_path || $user_database_path) {
  die "No reference genome was selected" 
}

# Allow over-ride of system-level database path with user
if ($user_database_path) {
  $database_path = $user_database_path;
  unless (`grep \\> $database_path`) {
      die "Error: $database_path  the user supplied file is not a FASTA file";
  }
  my $name = basename($database_path, qw/.fa .fas .fasta .fna/);
  print STDERR "hisat2-indexing $name\n";
  my $hisat2b = "hisat2-build";
  system $hisat2b . " $database_path $name";
  
  for my $query_file (@file_query) {
    # Grab any flags or options we don't recognize and pass them as plain text
    # Need to filter out options that are handled by the GetOptions call
    my @args_to_reject = qw(-xxxx);

    my $second_file = shift @file_query2 if @file_query2;

    my $HISAT_ARGS = join(" ", @ARGV);
    
    foreach my $a (@args_to_reject) {
        if ($HISAT_ARGS =~ /$a/) {
            report("Most TopHat arguments are legal for use with this script, but $a is not. Please omit it and submit again");
            exit 1;
        }
    }
    
    my $success = undef;
    system "mkdir output";

    my $app  = "/usr/bin/hisat2";
    my $format = $file_type;

    chomp(my $basename = `basename $query_file`);
    $basename =~ s/\.\S+$//;

    if ($format eq "PE") {
    	my $align_command = "$app $HISAT_ARGS -x $name -1 $query_file -2 $second_file -p $n_cores | samtools view -bS - > $query_file.bam";
    	report("Executing: $align_command\n");
    	system $align_command;
    	system "samtools sort $query_file.bam $basename.sorted";
    	system "mv $basename.sorted.bam  output";
    	system "samtools index output/$basename.sorted.bam";
    	system "rm -rf *bam";
        }
    elsif($format eq "SE"){
    	my $align_command = "$app $HISAT_ARGS -x $name -U $query_file -p $n_cores | samtools view -bS - > $query_file.bam";
    	report("Executing: $align_command\n");
    	system $align_command;
    	system "samtools sort $query_file.bam $basename.sorted";
    	system "mv $basename.sorted.bam  output";
    	system "samtools index output/$basename.sorted.bam";
    	system "rm -rf *bam";
        }
   }

}
system "rm -rf *.ht2";

sub report {
    print STDERR "$_[0]\n";
}
