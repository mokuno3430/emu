#!/usr/bin/perl

use warnings;
use strict;

use constant WINDOW => $ARGV[3];
use constant STEP => $ARGV[4];

&main;

sub main
{
    ( @ARGV != 5 ) and die "\nusage: $0 samtools-depth.out reference.fa read-depth_per_haploid window-size step-size > out.tsv\n\n";
    my ( %seq, %depth, %repeat );
    %seq = &input_fasta;
    &input_mpileup( \%seq, \%depth, \%repeat );
    &cal_copynumber( \%seq, \%depth, \%repeat );
}


sub input_fasta
{
    my $ID;
    my %seq;
    open IN, $ARGV[1];
    while( <IN> )
    {
	chomp;
	if( $_ =~ />(\S+)/ ){
	    $ID = $1;
	}else{
	    $seq{ $ID } .= $_;
	}
    }
    close IN;
    return %seq;
}


sub input_mpileup
{
    my $seq = shift;
    my $depth = shift;
    my $repeat = shift;
    my @buf;
    my $cov;
    open IN, $ARGV[0];
    while( <IN> )
    {
	chomp;
	@buf = split/\t+/;
	next if( ! exists $seq->{ $buf[0] } );
	$cov = 0;
	for( my $i = 2; $i <= $#buf; $i++ )
	{
	    $cov += $buf[$i];
	}
	if( $cov < $ARGV[2] * 4 ){ #count read depth at each position excluding repeat reagions
	    $depth->{$buf[0]}[int($buf[1]/STEP)] += $cov;
	}else{
	    $repeat->{$buf[0]}[int($buf[1]/STEP)]++;
	}
    }
    close IN;
}


sub cal_copynumber
{
    my $seq = shift;
    my $depth = shift;
    my $repeat = shift;
    my( $str, $read_depth, $repeat_region, $N_region, $valid_region, $ave );
    my $j = 0;
    foreach my $chr ( sort { length( $seq->{$b} ) <=> length( $seq->{$a} ) } keys %{$seq} )
    {
	last if( length( $seq->{$chr} ) < WINDOW );
	next if( ! exists $depth->{$chr} );
	for( my $i = 0; $i <= $#{$depth->{$chr}} - WINDOW/STEP; $i++ )
	{
	    $str = $i*STEP;
	    $read_depth = 0;
	    $repeat_region = 0;
	    $N_region = 0;
	    $valid_region = 0;
	    for( my $l = 0; $l < WINDOW/STEP; $l++ )
	    {
		$read_depth += $depth->{$chr}[$i+$l] if( $depth->{$chr}[$i+$l] );
		$repeat_region += $repeat->{$chr}[$i+$l] if( $repeat->{$chr}[$i+$l] );
	    }
	    $N_region++ while(  substr( $seq->{ $chr }, $str, WINDOW ) =~ m/N|n/g );
	    $valid_region = length( substr( $seq->{ $chr }, $str, WINDOW ) )- $N_region - $repeat_region;
	    if( $ARGV[3] / 10 < $valid_region ){
		$ave = $read_depth / ( $valid_region * $ARGV[2] );
	    }else{
		$ave = 0;
	    }
	    printf( "%s\t%d\t%d\t%d\t%.2f\n", $chr, length( $seq->{ $chr } ), $str, $j, $ave );
	    $j += STEP;
	}
    }
}
