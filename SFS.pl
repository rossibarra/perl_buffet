#!/usr/bin/perl

my $license="
 ########################################################################
#  Copyright (C) 2010 Jeffrey Ross-Ibarra <rossibarra@ucdavis.edu>       #
#  and Reed Cartwright <reed@scit.us>                                    #
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
#  Thanks are due to K. Thornton for insightful discussion and           #
#  A.J. Eckert for catching bugs in an earlier version.                  #
 ########################################################################
";

use strict;

my $version="2.2";

# calculates SFS from polytable files

# command line args
my %files = (); my $header=1; my $unround=0; my $condition=0; my $unfolded=1; my $verb_freq=0; my $num_bins=0;
#print "how: $#ARGV\n";
for my $arg (0..$#ARGV){ 
	if($ARGV[$arg] eq "-F"){ $arg+=1; 
		my $num=$ARGV[$arg]; 
		for my $numit (1..$num){ 
			#print "this: ",$ARGV[$arg+$numit],"\n";
			$files{$ARGV[$arg+$numit]} = 1 foreach(grep(-e, glob($ARGV[$arg+$numit])))
		}
	}
	if($ARGV[$arg] eq "-q"){ $header=0; }
	if($ARGV[$arg] eq "-f"){ $verb_freq=1; }
	if($ARGV[$arg] eq "-n"){ $arg+=1; $num_bins=$ARGV[$arg]; }
	if($ARGV[$arg] eq "--unrounded"){ $unround=1; }
	if($ARGV[$arg] eq "--folded"){ $unfolded=0; }
	if($ARGV[$arg] eq "--conditional"){ $condition=1; }
	if($ARGV[$arg] eq "--version"){ 
		die "SFS version $version\n$license\n";
	}
}	
my @files = keys(%files);
die "\n\nUsage:\tSFS.pl -F <number of files> <filenames> 

the <filenames> argument accepts the wildcard * character for multiple files

Optional command line arguments:

-q removes header (only if -f is not also chosen)
-n new sample size
-f prints out frequency counts of derived mutations at every site 
--condition conditions on observing a site as polymorphic in the new sample size
--unrounded returns unrounded frequencies
--folded calculates folded SFS
--version prints version and license information
\n\n" unless(@files );

foreach my $file (@files){

	#open file, get data
	open TABLE, "<$file";
	my @lines=<TABLE>;
	close TABLE;
	
	#get positions, sampsize, numsites ancestral state
	$lines[0]=~m/(\d+)\s+(\d+)\n/; my $sample_size = $1; my $num_sites = $2;
	next if $num_sites==0;
	chomp $lines[1];
	chomp $lines[2];
	
	#set num_bins=sample_size
	$num_bins = $sample_size if $num_bins == 0;
	
	my @positions; my @ancestral; 
	if($lines[1]=~m/^\D+\s+(\d+.*)/){ @positions=split(/\s+/,$1); }
	else{ @positions=split(/\s+/,$lines[1]); }
	if( $lines[2]=~m/^\S\t/){  $lines[2]="\t".$lines[2]; @ancestral = split(/\t/,$lines[2]);}
	else{ @ancestral = split(/\s+/,$lines[2]); }
	
	# get snp data for each individual
	my @snps; my @bob;
	for my $i ( 3..$sample_size+2 ){
		chomp $lines[$i];
		@bob=split(/\t/,$lines[$i]);
		push @snps, [ @bob ];
	}
	#calculate SFS
	my %SFS = &get_SFS( \@ancestral, \@snps, $num_sites, $sample_size, \@positions, $file );

	if ($num_bins == $sample_size){ 
		&printout( \%SFS, 1, $num_bins-1  ) if $condition;
		&printout( \%SFS, 0, $num_bins ) unless $condition;
	}
	else{
		my %newSFS=&prob( $num_sites, $sample_size, $num_bins, \%SFS );
		&printout( \%newSFS, 1, $num_bins-1  ) if $condition;
		&printout( \%newSFS, 0, $num_bins ) unless $condition;
	}
}	

sub prob {
	my( $num_snps, $sample_size, $num_bins, $sfs )=@_;
	my %SFS=%{$sfs};
	my %newSFS; for my $n (0..$num_bins){ $newSFS{$n}=0; } #initialize array
	my %SFSbins; #SFS
 	my %contotal; for my $j (1..$sample_size-1){ $contotal{$j} = $condition ? 0 : 1; } # total prob for conditional prob. calculation... start form 0 and add if conditional, set to 1 otherwise.
	
	#method from Nielsen R, Bustamante C, Clark AG, Glanowski S, Sackton TB, et al. 2005 
	#A Scan for Positively Selected Genes in the Genomes of Humans and Chimpanzees. 
	#PLoS Biol 3(6): e170. doi:10.1371/journal.pbio.0030170
	if( $num_bins < $sample_size ) {
  		for my $j (1..$sample_size-1){	
			for my $i (0..$num_bins){	
 				$SFSbins{$j}->[$i]=( ( &binom($sample_size - $j, $num_bins - $i) * &binom( $j, $i ) ) ) / &binom( $sample_size, $num_bins );	#eq. 1 from above paper.
 				if( $i > 0 && $i < $num_bins && $condition ){ $contotal{$j}+=$SFSbins{$j}->[$i] } # add up only nonfixed probs.
 			}
 		}
	}
	
	else {
 		for my $j (1..$sample_size-1){
			for my $i (0..$num_bins){			
				$SFSbins{$j}->[$i]=(&binom($sample_size,$j)*&binom($num_bins,$i)/&binom($num_bins+$sample_size,$j+$i))*($sample_size+1)/($num_bins+$sample_size+1);	
				if( $i > 0 && $i < $num_bins && $condition ){ $contotal{$j}+=$SFSbins{$j}->[$i] } # add up only nonfixed probs.
			}
		}
	}
	
	#divide by total prob (1 or sum of poly) and multiply by number
	for my $j (1..$sample_size-1){
		for my $i (0..$num_bins){	
			$SFSbins{$j}->[$i]=$SFSbins{$j}->[$i]/$contotal{$j};  # make probability relative
			my $multiplier = $SFS{$j} ? $SFS{$j} : 0; #number of SNPs at given freq.
			$newSFS{$i}+= $SFSbins{$j}->[$i]*$multiplier;
		}
	}
	return( %newSFS )
}
	
sub printout {
	my( $sfs, $start, $end  ) = @_;
	my %SFS=%{$sfs};
	for my $i ($start .. $end){
		print "$i\t" if $header;
	}
	print "\n" if $header;
	for my $i ($start .. $end){
		if( $SFS{$i} ){
			$SFS{$i}=int( $SFS{$i}*1000 )/1000 unless $unround;
			print "$SFS{$i}\t";
		}
		else{ print "0\t"; }
	}
	print "\n";
}	
	
sub get_SFS{
	my ($ank, $snp, $numsites, $samsize, $pos, $file  )=@_;
	my @ank=@{$ank};
	my @snp=@{$snp};
	my @errors;
	my %SFS=();
	my @positions=@{$pos};
	my $yime=localtime();
	if( $verb_freq == 1 ){ print "\nSFS version $version run on $yime\n$file\nobserved n=$samsize\tnew n: $num_bins\tS=$numsites\n"; }
	for my $site (1..$numsites){
		my $der = 0; my $anc = 0;
		if( $unfolded==0 ){ # calculate folded (unfolded is automatic from ancestral line)
			my $counter=0; 
			$counter++ while $snp[$counter][$site] eq "-" || $snp[$counter][$site] eq "?" || $snp[$counter][$site] eq "N";
			$ank[$site] = $snp[$counter][$site];
		}
		else{ # for unfolded
			next if $ank[$site] eq "?" || $ank[$site] eq "N" || $ank[$site] eq "-";
		}
		my $temp=0;
		for my $dude (0..$samsize-1){
			if($snp[$dude][$site] eq "N" || $snp[$dude][$site] eq "-" ){ next; }
			elsif($snp[$dude][$site] eq $ank[$site]){ # site is same as ancestor
				$anc++;				
			}
			elsif ( $unfolded==0){ $der++; }	# folded and sites != 'ancestor'
			elsif ( $der>0 && $snp[$dude][$site] eq $snp[$temp][$site] ){ $der++; $temp=$dude; } # unfolded and same as derived state
			elsif( $der>0 ){  next; } # unfolded, not ancestral, but 3rd nucleotide so ignored
			else{ $der++; $temp=$dude; } # if no derived yet, make first derived site
		}
		if( $verb_freq == 1){ 
			if( $unfolded == 1 ){ print "$positions[$site-1]\t$der\n"; }
			else{ $der < $samsize-$der ? print "$positions[$site-1]\t$der\n" : print "$positions[$site-1]\t",$samsize-$der,"\n"; }
		}
		if($der<=$anc && $unfolded == 0 ){ $SFS{$der}++; push(@errors,$site) if $der==0;}
		elsif( $unfolded == 0 ){ $SFS{$anc}++; push(@errors,$site) if $anc==0; }
		else{ $SFS{$der}++; push(@errors,$site) if $der==0; }
	}
	
	#prints header unless quieted, prints errors "?" if any exist
	if($header){ 
		print "\nSFS version $version run on $yime\n$file\nobserved n=$samsize\tnew n: $num_bins\tS=$numsites\n" unless $verb_freq;
		if( $SFS{0} ){ 
			print "?=$SFS{0}\n"; 
			print "sites with errors at: @errors\n";
		}
		else{ print "\n"; }
	} 
	return( %SFS );
}
	
#binomial coefficient
sub binom {
	my ($n, $k) = @_;
	if( $n < $k ){ return( 0 ); }
	else{	return( &fac($n)/( &fac($k)*&fac($n-$k) ) ); }
}

sub fac {
	my ($n) = @_;
	$n < 2 ? return 1 : return $n * fac($n-1);
}
