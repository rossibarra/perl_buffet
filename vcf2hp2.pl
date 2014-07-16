#! /usr/bin/env perl 
use strict;
use warnings;
use List::Util qw(sum);
use Getopt::Long;

my $tDP  = 10;     # default threshold for read-depth
my $tGQ  = 20;    # default threshold for sample genotype-quality
my $help = 0;

&GetOptions(
             'd:i'    => \$tDP,
             'q:i'    => \$tGQ,
             'h|help' => \$help,
);

if ($help){
    _help();
    exit;
}

my %IUPAC = (
              "0"  => "0",
              "N"  => "N",
              "++" => "+",
              "--" => "-",
              "+-" => "=",
              "-+" => "=",
              "AA" => "A",
              "CC" => "C",
              "GG" => "G",
              "TT" => "T",
              "AG" => "R",
              "CT" => "Y",
              "CG" => "S",
              "AT" => "W",
              "GT" => "K",
              "AC" => "M",
);

my @samples;

while (<STDIN>) {
    next unless ( /\#CHROM/ || !/\#/ );
    chomp;

    my @array = split /\t/;

    # samples ids
    if (/\#CHROM/) {
        @samples = @array[ 9 .. $#array ];

        #        for ($i =9; $i <= $#array; $i++){
        #            $samples{$i} = $array[$i];
        #        }
        #
		# print header
		print "rs#","\t",
			  "alleles","\t",
			  "chrom","\t",
			  "pos","\t",
			  "strand","\t",
			  "assembly#","\t",
			  "center","\t",
			  "protLSID","\t",
			  "assayLSID","\t",
			  "panelLSID","\t",
			  "QCcod","\t",
			  join ("\t",@samples),"\n";

        next;
    }

    my $chr = $array[0];
    my $pos = $array[1];
    my $id  = $array[2];

    my @alleleArray;
    push @alleleArray, $array[3];    # first allele in array is the reference
                                     # then comes the alternate alleles
    push @alleleArray, ( split /\,/, $array[4] );

    my $altqual    = $array[5];
    my $filter     = $array[6];
    my $infostring = $array[7];

    # is this position an indel?
    my $indel = ( $infostring =~ /INDEL/ ) ? 1 : 0;
    my @indelArray;
    if ($indel) {
        if ( length( $alleleArray[0] ) > length( $alleleArray[1] ) ) {
            push @indelArray, '+';
            push @indelArray, '-';
        }
        else {
            push @indelArray, '-';
            push @indelArray, '+';
        }
    }

    # format field
    my @format = split /\:/, $array[8];

    # search for DP,GT and GQ fields
    my ( $iGT, $iGQ, $iDP );
    for ( my $i = 0 ; $i <= $#format ; $i++ ) {
        if ( $format[$i] eq 'GT' ) {
            $iGT = $i;
        }
        elsif ( $format[$i] eq 'GQ' ) {
            $iGQ = $i;
        }
        elsif ( $format[$i] eq 'DP' ) {
            $iDP = $i;
        }
    }
    unless ( defined $iGT ) {
        print STDERR "GT field missing in format in line $chr:$pos\n";
        exit;
    }
    unless ( defined $iGQ ) {
        print STDERR "GQ field missing in format in line $chr:$pos\n";
        exit;
    }
    unless ( defined $iDP ) {
        print STDERR "DP field missing in format in line $chr:$pos\n";
        exit;
    }

    my @altquals = 0;
    my @gtquals  = 0;
    my @hp2Genotypes;

    for ( my $j = 9 ; $j <= $#array ; $j++ ) {
        my $sampStr = $array[$j];
        my @tmp     = split ":", $sampStr;
        my $sGQ     = $tmp[$iGQ];
        my $sDP     = $tmp[$iDP];

        my $call;

        #if no read depth, assign null call (N);
        if ( $sDP == 0 ) {
            $call = 'N';

            # if position doesn't meet thresholds, assign 'X'
        }
        elsif ( $sGQ < $tGQ || $sDP < $tDP ) {
            $call = 'X';
        }
        else {
            my ( $a1, $a2 ) = $tmp[$iGT] =~ /(\d)\/(\d)/;
            my $sGT;

            # if position is an INDEL
            if ($indel) {
                $sGT = $indelArray[$a1] . $indelArray[$a2];
            }
            else {

                #$sGT = $alleleArray[$a1] . $alleleArray[$a2];
                $sGT = join '', sort split '',
                  ( $alleleArray[$a1] . $alleleArray[$a2] );
            }
            $call = $IUPAC{$sGT};
            unless ($call) {
                print STDERR
"Can't translate genotype call ($sGT) for $chr:$pos $samples[$j]\n";
                exit;
            }

            if ( $sGT ne $alleleArray[0] || $sGT ne $indelArray[0] ) {
                push @altquals, $sGQ;
            }
            push @gtquals, $sGQ;
        }
        push @hp2Genotypes, $call;
    }

    my $aveAltQual = sprintf( "%.0f", sum(@altquals) / scalar(@altquals) );
    my $aveGTQual  = sprintf( "%.0f", sum(@altquals) / scalar(@gtquals) );

    print $id, "\t",
      join( "/", @alleleArray ), "\t",
      $chr, "\t",
      $pos, "\t",
      "+",        "\t",
      "AGPv2",    "\t",
      "MaizeDiv", "\t",
      "SBS",      "\t",
      "JeffRI",   "\t",
      "JeffRI",   "\t",
      $aveAltQual, "\t",
      join( "\t", @hp2Genotypes ), "\n";

}

sub _help {
    print STDERR<<HELP
USAGE: $0 -d $tDP -q $tGQ 

    [options]
    -d  => minimum read-depth in sample to make a call. Defaults to $tDP
           If minimum read-depth is not met, sample is assigned the 'N' genotype

    -q  => mininum genotype quality of sample to make a call. Defaults to $tGQ;
           If minimum quality is not met, sample is assigned the 'X' genotype

HELP
    
}
