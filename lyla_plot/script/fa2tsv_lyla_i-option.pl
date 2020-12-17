#!/usr/bin/perl

while( <> )	
{
    chomp;
    if( $_ =~ />(\S+)/ )
    {
	$ID = $1;
	push( @order, $ID );
    }else{
	$length{ $ID } += length( $_ );
    }
}

for( my $i = 0; $i <= $#order; $i++ )
{
    printf( "%s\t0\t%d\t+\t%s\n", $order[$i], $length{ $order[$i] }, $order[$i] );
}
