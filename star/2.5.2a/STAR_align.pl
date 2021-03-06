#!/usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);


my (@file_query, $database_path, $user_database_path, $annotation_path,
$user_annotation_path, @file_query2, $file_type);


GetOptions( "file_query=s"      => \@file_query,
            "file_query2=s"     => \@file_query2,
            "user_database=s"   => \$user_database_path,
            "user_annotation=s" => \$user_annotation_path,
	          "file_type=s"       => \$file_type,
            );

# sanity check for input data
if (@file_query2) {
    @file_query && @file_query2 || die "Error: At least one file for each paired-end is required\n";
    @file_query == @file_query2 || die "Error: Unequal number of files for paired ends\n";
}

if (!($user_database_path)) {
    die "No reference genome was supplied\n";
}
if (!($user_annotation_path)) {
    die "No reference annotation was supplied\n";
}
if (@file_query < 1) {
    die "No FASTQ files were supplied\n";
}
# Sanity check for input ref. genome and annotation
unless ($user_database_path) {
  die "No reference genome was selected"
}

unless ($user_annotation_path) {
die "No reference genome was selected"
}

# Allow over-ride of system-level database path with user
if ($user_database_path) {
  $database_path = $user_database_path;
  unless (`grep \\> $database_path`) {
      die "Error: $database_path  the user supplied file is not a FASTA file";
  }
  my $name = basename($database_path, qw/.fa .fas .fasta .fna/);
  print STDERR "STAR-indexing $name\n";
  system "mkdir index";
  my $STARp = "STAR";
  system "$STARp --runThreadN 4  --runMode genomeGenerate  --genomeDir index --genomeFastaFiles $database_path --sjdbGTFfile $user_annotation_path";
}

for my $query_file (@file_query) {
    # Grab any flags or options we don't recognize and pass them as plain text
    # Need to filter out options that are handled by the GetOptions call
    my @args_to_reject = qw(-xxxx);


    my $second_file = shift @file_query2 if @file_query2;

    my $STAR_ARGS = join(" ", @ARGV);
    foreach my $a (@args_to_reject) {
        if ($STAR_ARGS =~ /$a/) {
            report("Most STAR arguments are legal for use with this script, but $a is not. Please omit it and submit again");
            exit 1;
        }
    }
    my $app  = "STAR";

    my $format = $file_type;
    
    chomp(my $basename = `basename $query_file`);
    $basename =~ s/\.\S+$//;
	 if ($format eq 'PE') {
    my $align_command = "$app $STAR_ARGS --runThreadN 4 --genomeDir index --outFileNamePrefix $basename --readFilesIn $query_file $second_file --readFilesCommand gunzip -c --sjdbGTFfile $user_annotation_path";
    
    report("Executing: $align_command\n");
    system $align_command;
    my $move_files = "mkdir output;mkdir output\$basename;mv Log* $basename* output\$basename";
    system $move_files;
	}
    elsif($format eq 'SE'){
	 my $align_command = "$app $STAR_ARGS --runThreadN 4 --genomeDir index --outFileNamePrefix $basename --readFilesIn $query_file --readFilesCommand gunzip -c --sjdbGTFfile $user_annotation_path";

    report("Executing: $align_command\n");
    system $align_command;
    my $move_files = "mkdir output;mkdir output\$basename;mv Log* $basename* output\$basename";
    system $move_files;

	}
}

sub report {
    print STDERR "$_[0]\n";
}
