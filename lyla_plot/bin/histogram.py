import matplotlib.pyplot as plt
import sys
import csv
import common
import os.path


def plot_histogram( seqs, ax, size, fn ):
    min_y = 0
    max_y = 4
    color = 'tomato'
    msize=3
    
    with open( fn , 'r' ) as file:
        x_list = []
        y_list = []
        for line in file:
            buf = line.rstrip( '\n' ).split( '\t' )
            i = common.detect_index( buf[0], int( buf[2] ), int( buf[2] ), seqs )
            if( i == -1 ):
                continue
            if( float( buf[4] ) < min_y or max_y < float( buf[4] )):
                continue
            x = seqs[i][buf[0]].convert_position2xcoord( int(buf[2] ))
            y0 = seqs[i][buf[0]].origin_y + seqs[i][buf[0]].height + size.margin_bw_scaffold_alignment
            y1 = seqs[i][buf[0]].origin_y + size.histograms[i] - size.margin_bw_scaffold_alignment
            y = ( (y1 - y0 ) * (float( buf[4] ) - min_y ) )/(max_y - min_y) + y0
            x_list.append( x )
            y_list.append( y )
    if( len( x_list ) != 0 ):
        ax.scatter( x_list, y_list, color=color, marker='.', s=msize, zorder=2, lw=0 )
            
            
def plot_background( seqs, ax, size ):
    color = common.Color( 'gainsboro', 0.5 )    
    for i in range( len( seqs )):
        if( size.histograms[i] == 0 ):
            continue
        for scaf in seqs[i]:
            x_start = seqs[i][scaf].origin_x
            x_end = seqs[i][scaf].origin_x + seqs[i][scaf].length
            y_bottom = seqs[i][scaf].origin_y + size.margin_bw_scaffold_alignment + seqs[i][scaf].height
            y_top = seqs[i][scaf].origin_y + size.histogram_height - size.margin_bw_scaffold_alignment
            x = [ x_start, x_end, x_end, x_start ]
            y = [ y_bottom, y_bottom, y_top, y_top ]
            seqs[i][scaf].plot_line_on_histogram( ax, size.histograms[i] )
            if seqs[i][scaf].name == 'BLANK':
                continue
            ax.fill( x, y, color=color.color, lw=0, alpha=color.alpha, zorder=1 )
        
