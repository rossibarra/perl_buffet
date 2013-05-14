#!/usr/bin/perl

use strict;
use Getopt::Std;
# get options and input files
my %options=();
getopts("i:vcmgntqsp",\%options) or help();
change() if $options{c};

# gets infiles
my %files = ();
$files{$options{i}} = 1 foreach(grep(-e, glob($options{i}))); 
my @files = keys(%files);


foreach my $file (@files){

	open FILE, "<$file";

	my %seqhash=(); my $alignment; my $modified=0; my $seq;
	while( <FILE>){
		if( $modified ){
 			if($_=~m/\<alignment\>(.*)<\/alignment\>/){
 				$seq = $1;
			}
			elsif( $_=~m/\<name\>(.*)<\/name\>/ ){
				my $ind = $1; $seqhash{$ind}=$seq;
			}
			elsif( $_=~m/\<\/sequences\>/ ){
				$modified=0;
			}
 		}
 		else{
			if( $_=~m/\<name\>(.*)\<\/name\>/ ){ 
				if( $alignment ){
					open OUT, ">$alignment.fasta";
					foreach my $dude ( keys(%seqhash) ){
						print OUT ">$dude\n$seqhash{$dude}\n";
					}
					close OUT;
				}
				%seqhash=();
				$alignment=$1;
				$modified=1 if $_=~m/\(modified\)/;
			}
			next if $_=~m/\<|\*|^\s/;
			my $data=$_; chomp $data;
			my @temp=split(/\s+/,$data);
			my $ind=$temp[0];
			$seq=$temp[1];
			if( $seqhash{$ind} ){ $seqhash{$ind}=$seqhash{$ind} . $seq ; }
			else{ $seqhash{$ind} = $seq; }
		}
	}
		
	open OUT, ">$alignment.fasta";
	foreach my $dude ( keys(%seqhash) ){
		print OUT ">$dude\n$seqhash{$dude}\n";
	}
	close OUT;
	
}
