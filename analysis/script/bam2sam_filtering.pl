#!/usr/bin/perl

use warnings;
use strict;

( @ARGV != 3 ) and die "\nusage: $0 single-mode_raw.bam min_identity(%) min_alignment_cover_rate(%) > out.sam\n\n";


my ( @buf, @match );
my $len;
my $SAMTOOLS = 'samtools';

open IN, "$SAMTOOLS view -H $ARGV[0] |";
while( <IN> )
{
    print $_;
}
close IN;

open IN, "$SAMTOOLS view $ARGV[0] |";
while ( <IN> )
{	
    chomp;
    @buf = split/\t/;
    next if( $buf[1] > 16 );
    @match = ( $buf[5] =~ m/(\d+)M/g );
    $len = 0;
    for( my $i = 0; $i <= $#match; $i++ )
    {
	$len += $match[$i];
    }
    next if( $len / length( $buf[9] ) *100 < $ARGV[2] ); #alignment cover rate
    next if( $_ !~ /NM:i:(\d+)/ );	
    next if( ( $len - $1 )/$len * 100 < $ARGV[1] ); #identity
    print $_,"\n";
}
close IN;


