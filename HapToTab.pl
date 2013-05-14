#!/usr/bin/perl
use strict;
use warnings;

if( $ARGV[0] eq "--help" || $#ARGV<1 ){ die "
 ########################################################################
#  Copyright (C) 2011 Jeffrey Ross-Ibarra <rossibarra\@ucdavis.edu>      #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  (at your option) any later version.                                   #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#  GNU General Public License <http://www.gnu.org/licenses/>             #
#  for more details.                                                     #
#                                                                        #
 ########################################################################

	HapToTab version 2.1 can be run as follows:
	
	HapToTab.pl <INFILE> <--snp or --bp> <INT> <OPTIONS>
	
	Where <INFILE> is a file in hapmap format.  
	HapToTab needs to know the window size either in SNPs or bp (e.g. --snp 100 or --bp 10000).  
	A negative number for SNPs will run HapToTab using only one window for the whole file.  
	The following options are available:
	
	--g <NAME> name the group you are running haptotab on.  default is \"all\"
	--maxdip format the output files for Hudson's maxdip program
	--maxhap format the output files for Hudson's maxhap program
	--nohets remove all heterozygotes (code as missing data)
	--out file has an outgroup (assumes the last column in the hapmap file)
	--freq <INT> filters out all SNPs with a frequency less than INT in the data\n\n";
}
open FILE, "<$ARGV[0]" or die "cannot open file '$ARGV[0]'!\n\n";

#define stuff
my $teomz="all";
my $window_by_snp=0;
my $maxdip=0;
my $snp_window_size;
my $bp_window_size;
my $outgroup=0;
my $splithets=1;
my $freqfilter=0;
my $filter_freq=999;
my $hapster=0;

for my $arg (1..$#ARGV){
	if( $ARGV[$arg] eq "--snp" ){ $window_by_snp=1; $snp_window_size=$ARGV[$arg+1]; }
	if( $ARGV[$arg] eq "--bp" ){ $window_by_snp=0; $bp_window_size=$ARGV[$arg+1]; }
	if( $ARGV[$arg] eq "--g" ){ $teomz=$ARGV[$arg+1];; }
	if( $ARGV[$arg] eq "--maxdip" ){ $maxdip=1; }
	if( $ARGV[$arg] eq "--maxhap" ){ $splithets=0; $maxdip=1; $hapster=1; }
	if( $ARGV[$arg] eq "--nohets" ){ $splithets=0; }
	if( $ARGV[$arg] eq "--out" ){ $outgroup=1; }
	if( $ARGV[$arg] eq "--freq" ){ $freqfilter=1; $filter_freq=$ARGV[$arg+1]; }
}

#grab header, ditch first 11 columns to get only individual line names
#assumes infile has been tweaked with cut to include only correct lines
#assumes outgroup, if present, is LAST column in file
my ($header, @lines)=<FILE>; close FILE;
chomp $header;
my @inds=split(/\t/,$header);
for my $i (1..11){
	shift(@inds);
}
pop(@inds) if $outgroup;

die "You have too many individuals.  This will cause some errors.  Suck it up and fix this code.\n\n" if $#inds>997;

my $num_inds=$#inds+1;
my $snp_counter=0; #snp iterator
my $bp_start=0; #window bp iterator
my $bp_counter=0; #window bp length
my %window=();
my $snp;
my $csome;
my $pos;
my @alleles;
my $win=1;
my $pos_old=0;

#run through lines
for my $L ( 0..$#lines ){	

	my $line=$lines[$L];
	next if $line=~m/\*/ || $line=~m/\-\/\+/ || $line=~m/\+\/\-/; # all indels appear to have a "*" or "-/+" combo in the allele, so regex 'em out.	
	chomp $line;  #you know it
		
	#regex to grab alleles 
	$line=~m/(\S+)\t\S+\t(\d+)\t(\d+)\t\S\t.*\t\d+\t(\w+\t.*)/;
	$snp=$1;
	$csome=$2;
	$pos=$3;		
	@alleles=split(/\t/,$4);
	
	#skips if a site is called twice
	if( $L<$#lines ){ 
		my $next_line=$lines[$L+1];
		$next_line=~m/(\S+)\t\S+\t(\d+)\t(\d+)\t\S\t.*\t\d+\t(\w+\t.*)/;	
		my $next_pos=$3;
		my $last_line=$lines[$L-1];
		$last_line=~m/(\S+)\t\S+\t(\d+)\t(\d+)\t\S\t.*\t\d+\t(\w+\t.*)/;	
		my $last_pos=$3;
		next if $pos == $next_pos || $pos == $last_pos;
	}
	
	#iterate stuff
	$bp_counter = $pos-$bp_start;
	$snp_counter++;
		
	#reached snp_window_size snps
	if( $window_by_snp == 1 && $snp_counter == $snp_window_size ){ 
		#since ended on this line, add this line, then print		
		$window{$pos}=[@alleles];
		&poly_make(\%window,$snp_counter, $bp_start+1, $pos, $snp_window_size );
			
		#clean shit up
		$bp_start=$pos; $bp_counter=0;  $snp_counter=0; %window=(); $pos_old=$pos;
	}
	#THIS line pushes over bp_window_size bp
	elsif( $window_by_snp == 0 && $bp_counter > $bp_window_size ){ 
		#since window ended on last line, print, then add this line
		&poly_make(\%window, $snp_counter-1, $bp_start+1, $bp_start+$bp_window_size, $win );
		$win++;
		#clean shit up
		while( $bp_start+$bp_window_size < $pos ){
			$bp_start=$bp_start+$bp_window_size; 
		}
		$bp_counter=0; $snp_counter=1; %window=();
		
		#add in this line after cleanup
		$window{$pos}=[@alleles];
		$pos_old=$pos;
	}
	#if window not ended, add stuff;
	else{ 
		$window{$pos}=[@alleles]; 
		$pos_old=$pos; 
	}
}
#push out last window only if window_size is negative (all snps are desired), otherwise drop these SNPs since window will not be comparable.
&poly_make(\%window,$snp_counter, $bp_start+1, $pos, $win ) if $window_by_snp==1 && $snp_window_size < 0;
&poly_make(\%window, $snp_counter, $bp_start+1,  $bp_start+$bp_window_size, $win ) if $window_by_snp==0;

#make a polytable
sub poly_make
{
	my ($window_ref, $num_snps, $start, $end, $win)=@_;
	my %window= %{ $window_ref };
	my @positions = sort {$a <=> $b} keys(%window);
	my $poly = $window_by_snp == 1 ? "poly.snp" : "poly.bp";
 
 	my $file="$poly.$teomz.$csome.$start.$end.$win";
 	open OUT, ">$file";
 
 	#get number of dudes for maxdip, normal
 	my $howmany = $maxdip == 1 ? $num_inds : 2*$num_inds;	
 
	#calculate MAF for each locus
	my @freqvector=();
	if( $freqfilter ){
		foreach my $p (@positions){
			my %allelehash=();
			my $denominator = $num_inds;  ## to account for things that are fixed
			foreach my $i ( 0..$#inds ){
				my $allele=$window{$p}->[$i];
				my $gamete=&disambiguity( $allele );
				if( $gamete=~m/N/ ){ $denominator--; next; }
				print STDERR $gamete, "\t", substr($gamete,0,1), "\n";
				if( $allelehash{substr($gamete,0,1)} ){ $allelehash{substr($gamete,0,1)}++ }
				else{	$allelehash{substr($gamete,0,1)}=1; }
				if( $splithets ){
					if( $allelehash{substr($gamete,1,1)} ){ $allelehash{substr($gamete,1,1)}++ }
					else{	$allelehash{substr($gamete,1,1)}=1; }
				}
			}
			my $min=$num_inds;
			foreach(keys(%allelehash)){ 
				$min=$allelehash{$_} if $allelehash{$_} < $min;
			#	print STDERR $_, "\t", $allelehash{$_}, "\t", $denominator, "\t", $min, "\n";
				$min=0 if $allelehash{$_} == $denominator; #added in case "SNP" is all N's and one nucleotide
			#	print STDERR $_, "\t", $allelehash{$_}, "\t", $denominator, "\t", $min, "\n";
			}
			$freqvector[$p]=$min;
			$num_snps-- if $freqvector[$p]<=$filter_freq;
		}
	}

 	#modified to remove SNPs below frequency
 	print OUT "$howmany\t$num_snps\n\t";
 	foreach my $p (@positions){ print OUT "$p\t" unless $freqfilter && $freqvector[$p]<=$filter_freq; }
 	print OUT "\nANC\t";
 	if( $outgroup ){
		foreach my $p ( @positions ){
			next if $freqvector[$p]<=$filter_freq; 
			my $allele=$window{$p}->[$#inds+1];
						
			$allele="N" if $allele=~m/[K|M|R|Y|S|W]/; 
			print OUT "$allele\t";
		}
	}
	else{
		foreach my $p (@positions){ print OUT "?\t" unless $freqfilter && $freqvector[$p]<=$filter_freq; }
	}
	print OUT "\n";
	
	#print to file
	for my $counter (0..$#inds){
		#take allele and make it diploid
		for my $copy (0..1){
			next if $copy == 1 && $hapster;
			print OUT "$inds[$counter]\t";
			foreach my $p ( @positions ){
				my $allele=$window{$p}->[$counter];
				my $gamete=&disambiguity( $allele );
			#	print STDERR "$freqvector[$p]\t$filter_freq\n";
				next if $freqfilter && $freqvector[$p]<=$filter_freq;
				print OUT substr($gamete,$copy,1),"\t";
			}
			print OUT "\n";
		}
	}
 	close OUT;
}
	
sub disambiguity
{
	my $allele=$_[0];

	#unambiguize them damn codes
	if($splithets){ $allele=~s/K/GT/; $allele=~s/M/AC/; $allele=~s/R/AG/; $allele=~s/Y/CT/; $allele=~s/S/CG/; $allele=~s/W/AT/; }
	else{ $allele=~s/K/N/; $allele=~s/M/N/; $allele=~s/R/N/; $allele=~s/Y/N/; $allele=~s/S/N/; $allele=~s/W/N/; }
		
	#duplicate if only one nuke
	length($allele) <2 && !$hapster ? return( $allele . $allele) : return($allele);
}
