 ##########################################################################
#																									#	
#  Copyright (C) 2008 Jeffrey Ross-Ibarra <rossibarra@gmail.com>     	  	#
#                                                                         	#
#  This program is free software: you can redistribute it and/or modify   	#
#  it under the terms of the GNU General Public License as published by   	#
#  the Free Software Foundation, either version 3 of the License, or      	#
#  (at your option) any later version.                                    	#
#                                                                         	#
#  This program is distributed in the hope that it will be useful,        	#
#  but WITHOUT ANY WARRANTY; without even the implied warranty of         	#
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          	#
#  GNU General Public License <http://www.gnu.org/licenses/>					#
# 	for more details.																			#
#																									#
 ##########################################################################

#!/usr/bin/perl
package main;

	use lib '/Users/jri/bin/bioperl-1.4/';  
	use strict;
	use Getopt::Std;
	use Bio::DB::GenBank;
	
	my $gb = Bio::DB::GenBank->new();
	my $version="0.1";
	
	# get options and input files
	my %options=();
	getopts("I:VOCNMSG",\%options) or help();
	
	# gets infiles
	my %files = ();
	$files{$options{I}} = 1 foreach(grep(-e, glob($options{I}))); 
	my @files = keys(%files);
	&Help unless(@files);
	
	#do stuff to each file
	foreach my $file (@files){
	
		open FILE, "<$file";	
		open TEMP, ">$file.one";
		my $name="";
		my $seq="";

		while(<FILE>){
			if( $_=~m/\>/ ){
				if( length($seq)>0 ){
					print OUT "$name\n$seq\n";
				}	
				$name=$_; 
				chomp $name;
				$seq="";
			}
			else{
				chomp $_;
				$seq.=$_;
			}
		}
		print TEMP "$name\n$seq\n";
		close TEMP;
		close FILE;
		$file="$file.one";

		#new fastA object
		my $bob = new FastA;
		
		#verbose output
		&Verb(\%options, $file) if ($options{V});
		
		#changelog
		&Change() if $options{C};
		
		# reads in a fastafile with EF names	
		$bob->read($file);
		
		#turns end gaps into N
		$bob=&Ngaps($bob) if ($options{N});
		
		#turns weird characters into N
		$bob=&Maken($bob) if ($options{M});
		
		#makes genbank accesion numbers into species
		$bob=&Genbank($bob) if ($options{G});
		
		#prints to outfile or screen
		my $outfile = $file.".out";
		if ($options{O}){ &Out($bob,$outfile); }
		elsif( !$options{G} ){ &screen($bob); }
		
		#splits groups based on individuals ID's and write to separate files
		&Split($bob,$outfile) if $options{S};
	}
	
	#SUBS##
	
	#prints version and options used
	sub Verb{
			my %options=%{@_[0]};
			my $file=@_[1];
			my $yime=localtime();

			print "\n\nfastU (Fasta Utilies) version $version output for file $file on $yime\nOptions used: ";
			foreach my $used_opt ( keys(%options) ){ 
				print " -", $used_opt unless $used_opt=~m/i/;  
			}
			print "\n\n";	
	}
	
	#options print to stdout
	sub Help{
		print "\n********************************************\n*";
		print "\n* fastU $version Copyright 2008 Jeffrey Ross-Ibarra <rossibarra\@gmail.com>\n*";
		print "\n*\tfastU performs a number of fasta-handling utilities\n*";
		print "\n*\tRequired command line parameters:\n*\n*";
		print "\t\t-I 'filename' : input name of fasta file(s) (accepts * wildcard).\n*";
		print "\n*\tOptional command line parameters:\n*\n*";
		print "\t\t-V : for verbose output\n*";
		print "\t\t-N : turn end gaps into N's\n*";
		print "\t\t-M : turn all characters except A,C,G,T,_ into N's\n*";
		print "\t\t-O : outputs to file 'filename'.out instead of stdout\n*";
		print "\t\t-C : changelog since version 1.0\n*";
		print "\t\t-S : splits alignment into separate files for groups of sequence ID's\n*";
		print "\t\t-G : replaces Genbank ID's with species name (requires Bioperl instalation)\n*";
		print "\n********************************************\n";
		die "\n";
	}
	
	#changes since verion 0.1
	sub Change{
		print "\n********************************************\n*";
		print "\n* fastU $version Copyright 2008 Jeffrey Ross-Ibarra <rossibarra\@gmail.com>\n*";
		print "\n*\tChanges since v 0.1:\n*";
		print "\n********************************************\n";
		die "\n";
	}
	
	# turns end gaps into N's
	sub Ngaps{
		my $sue=@_[0];
		my @inds = keys%$sue;
		#print "@inds";
		foreach (@inds) {
			my $seq = $sue->get($_) or print "Cannot find sequence $_ in subroutine Ngaps\n";
			# beginning gap
			if( $seq=~m/(^\-+)/){
				my $gap=$1;
				my $gapN=$gap;  $gapN=~s/\-/N/g;
				$seq=~s/^$gap/$gapN/;
				$sue->set($_,$seq);
			}
			#end gap
			if( $seq=~m/(\-+$)/){
				my $gap=$1;
				my $gapN=$gap;  $gapN=~s/\-/N/g;
				$seq=~s/$gap$/$gapN/;
				$sue->set($_,$seq);
			}	
		}
		return($sue);
	}
	
	# turns all non ACGT- into N's
	sub Maken{
		my $sue=@_[0];
		my @inds = keys%$sue;
		foreach my $dude (@inds) {
			my $seq = $sue->get($dude) or print "Cannot find sequence $dude in subroutine Maken\n";	
			$seq=~s/[^A|C|T|G|N|\-]/N/g;
			$sue->set($dude,$seq);
		}
		return($sue);
	}
	
	# gets genbank species name if fasta file has accession numbers
	sub Genbank{
		my $sue=@_[0];
		my @inds = keys%$sue;
		my $ted = new FastA;
		foreach my $dude (@inds) {
			# print "dude: $dude\n";
			my $seq = $sue->get($dude) or print "Cannot find sequence $dude in subroutine Genbank\n";	
			my $gbseq=$gb->get_Seq_by_acc($dude) or print "$dude is not a genbank accession #!\n";
			my $species=$gbseq->species->binomial . " " . $gbseq->species->sub_species . " (".$dude.")";
			$ted->{$species}=$seq;
		}
		return($ted);
	}
	
	#prints fasta alignment to outfile
	sub Out{
		my $sue=@_[0];
		my $outfasta=@_[1];
		
		open OUT, ">$outfasta";
		my @inds = keys%$sue;
		foreach my $dude (@inds) {
			my $seq = $sue->get($dude) or print "Cannot find sequence $dude in sub Out\n";
			print OUT ">$dude\n$seq\n";		
		}
		close OUT;
	}
	
	sub Split{
		my $sue=@_[0];
		my $tempout=@_[1];
		my @inds = keys%$sue;
		
		foreach my $id (@inds){
			my $seq = $sue->get($id) or print "Cannot find sequence $id in sub Split\n";
			$id=~s/\s//g;
			$id=~m/(\w+)\(.*\)/;
			$id=$1;
			my $tempfile=$tempout.$id; 
			open TEMP, ">>$tempfile";
			print TEMP ">$id\n$seq\n";		
			close TEMP;
		}
	}
	
	sub screen{
		my $sue=@_[0];
	
		my @inds = keys%$sue;
		foreach my $dude (@inds) {
			my $seq = $sue->get($dude) or print "Cannot find sequence $dude in sub screen\n";
			print ">$dude\n$seq\n";		
		}
	}

package FastA;

	sub new {
		my ($class) = @_;
		my ($self) = {};
		bless ($self, $class);
		return $self;
	}
	
	sub read {
		my ($self, $file) = @_;
		my ($line, $ID);
	
		open (FASTAFILE, $file) or return 0;
	
		while ($line = <FASTAFILE>) {
	
			# sequence name
			if ($line =~ /^\>(.*)\n/) {
				$ID = $1;
			}
			# or sequence
			else {
				chomp $line;
				$self->{$ID} = $line;
			}
		}
	
	close FASTAFILE;
	return 1;
	}
	
	sub get {
		my ($self, $key) = @_;
	
	return $self->{$key};
	}
	 
	sub set {
		my ($self, $key, $value) = @_;
		$self->{$key} = $value;
	}
	
	1;
