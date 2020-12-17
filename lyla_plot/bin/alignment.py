import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as clr
import sys
import csv
import common
import os.path

class Colormap:
    pallet = {
        'red':((0.00, 0.10, 0.10),
               (0.25, 0.10, 0.10),
               (0.40, 0.30, 0.30),
               (0.60, 1.00, 1.00),
               (0.80, 0.90, 0.90),
               (1.00, 0.70, 0.70)),
        'green':((0.00, 0.10, 0.10),
                 (0.25, 0.60, 0.60),
                 (0.40, 0.80, 0.80),
                 (0.60, 0.75, 0.75),
                 (0.80, 0.30, 0.30),
                 (1.00, 0.15, 0.15)),
        'blue':((0.00, 0.40, 0.40),
                (0.25, 1.00, 1.00),
                (0.40, 0.25, 0.25),
                (0.60, 0.00, 0.00),
                (0.80, 0.05, 0.05),
                (1.00, 0.20, 0.20))
    }
    mycm = clr.LinearSegmentedColormap('original', pallet )

    cmaps = [ cm.binary, cm.bone_r, cm.inferno_r, cm.hot_r, cm.YlGnBu, mycm ]
    cmap_list = [ 'binary', 'bone_r', 'inferno_r', 'hot_r', 'YlGnBu', 'original' ]

    def __init__( self, min_identity, max_identity, cm, alpha=0.5 ):
        self.min_identity = min_identity
        self.max_identity = max_identity
        self.alpha = alpha
        self.cm = cm

    def convert_identity2color( self, identity ):
        color_value=( identity - self.min_identity )/( self.max_identity - self.min_identity )
        align_color=Colormap.cmaps[ self.cm ]( color_value )
        return align_color

    def output_parameters( self ):
        print( '##Colormap paramenters:' )
        print( '  min_identity: %d' % ( self.min_identity))
        print( '  max_identity: %d' % ( self.max_identity ))
        print( '  alpha: %.2f' % ( self.alpha ))
        print( '  colormap: %d (%s)' % ( self.cm, Colormap.cmap_list[self.cm]  ))
        print( '' )

    
class Colorbox:
    def __init__( self, size ):
        self.height = size.bottom_margin * 0.2
        self.width = size.xlim_max * 0.2
        self.origin_x = size.xlim_max * 0.55
        self.origin_y = float( ( size.bottom_margin - self.height ) /2 )
        
    def plot( self, ax, heatmap ):
        cell_width = self.width / ( heatmap.max_identity - heatmap.min_identity )
        text_legend = common.Text( 'identity (%)', 8, self.origin_x - cell_width * 2, self.origin_y + self.height / 2,  'right', 'center' )
        text_legend.output( ax )
        cell_originx = self.origin_x
        for i in range( heatmap.min_identity, heatmap.max_identity + 1 ):
            align_color=heatmap.convert_identity2color( i )
            colorbox_x = [ cell_originx, cell_originx + cell_width, cell_originx + cell_width, cell_originx ]
            colorbox_y = [ self.origin_y, self.origin_y, self.origin_y + self.height, self.origin_y + self.height ]
            cell_originx += cell_width
            ax.fill( colorbox_x, colorbox_y, color=align_color, alpha=heatmap.alpha, linewidth=0 )
            if( i % 10 == 0 ):
                num_legend = common.Text( i, 6, cell_originx, self.origin_y - 0.06, 'center', 'top' )
                num_legend.output( ax )

    def output_parameters( self ):
        print( '##Colorbox paramenters:' )
        print( '  width: %.2f' % ( self.width ))
        print( '  height: %.2f' % ( self.height ))
        print( '  origin_x %.2f' % ( self.origin_x ))
        print( '  origin_y: %.2f' % ( self.origin_y ))
        print( '' )

                
def convert_position2coord( A, B, A_start, A_end, B_start, B_end, margin, A_h_height, B_h_height ):
    x1 = A.convert_position2xcoord( A_start )
    x2 = A.convert_position2xcoord( A_end )
    x3 = B.convert_position2xcoord( B_end )
    x4 = B.convert_position2xcoord( B_start )    
    x = [ x1, x2, x3, x4 ]
    
    if( A.origin_y < B.origin_y ):
        y1 = A.convert_position2ycoord( margin + A.height + A_h_height )
        y2 = B.convert_position2ycoord( margin * -1 )
    else:
        y1 = A.convert_position2ycoord( margin * -1 )
        y2 = B.convert_position2ycoord( margin + A.height + B_h_height )
    y = [ y1, y1, y2, y2 ]

    return x, y


def count_alignment_files( args ):
    valid_files = 0
    input_formats = [ args.alignment, args.blastn, args.lastz, args.mummer ]
    for files in input_formats:
        if files is None:
            continue
        for fn in files:
            if not os.path.isfile( fn ):
                continue
            valid_files += 1
    return valid_files


def set_min_identity( args ):
    list_min_identity = []
    input_formats = [ args.alignment, args.blastn, args.lastz, args.mummer ]
    func_set_min_identity = [ cal_min_identity4original, cal_min_identity4blastn, cal_min_identity4lastz, cal_min_identity4mummer ]
    if( args.min_identity != -1 ):
        return args.min_identity
    for files, func_min in zip( input_formats, func_set_min_identity ):
        if files is None:
            continue
        for fn in files:
            if not os.path.isfile( fn ):
                print( 'WARNING: %s is not found\n' % fn )
                continue
            list_min_identity.append( func_min( fn ) )
    if( len( list_min_identity ) == 0 ):
        return 0
    return min( list_min_identity )


def plot_alignment4original( seqs, ax, heatmap, size, fn ):
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( '\t' )
            i = common.detect_index( buf[0], int( buf[1] ), int( buf[2] ), seqs )
            if( buf[4] < buf[5] ):
                j = common.detect_index( buf[3], int( buf[4] ), int( buf[5] ), seqs )
            else:
                j = common.detect_index( buf[3], int( buf[5] ), int( buf[4] ), seqs )
            color = heatmap.convert_identity2color( float( buf[6] ))

            if( i == -1 or j == -1 ):
                continue

            x, y = convert_position2coord( seqs[i][buf[0]], seqs[j][buf[3]], int( buf[1] ), int( buf[2] ), int( buf[4] ), int( buf[5] ), size.margin_bw_scaffold_alignment, size.histograms[i], size.histograms[j] )
            ax.fill( x, y, color=color, alpha=heatmap.alpha, lw=0 )

            
def plot_alignment4blastn( seqs, ax, heatmap, size, fn ):
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( '\t' )
            i = common.detect_index( buf[0], int( buf[6] ), int( buf[7] ), seqs )
            if( buf[8] < buf[9] ):
                j = common.detect_index( buf[1], int( buf[8] ), int( buf[9] ), seqs )
            else:
                j = common.detect_index( buf[1], int( buf[9] ), int( buf[8] ), seqs )
            color = heatmap.convert_identity2color( float( buf[2] ))

            if( i == -1 or j == -1 ):
                continue

            x, y = convert_position2coord( seqs[i][buf[0]], seqs[j][buf[1]], int( buf[6] ), int( buf[7] ), int( buf[8] ), int( buf[9] ), size.margin_bw_scaffold_alignment, size.histograms[i], size.histograms[j] )
            ax.fill( x, y, color=color, alpha=heatmap.alpha, lw=0 )

    
def plot_alignment4lastz( seqs, ax, heatmap, size, fn ):
    with open( fn , 'r' ) as file:
        for line in file:
            if( line[0:1] == "#" ):
                continue
            buf = line.rstrip( '\n' ).split( '\t' )
            i = common.detect_index( buf[1], int( buf[4] ), int( buf[5] ), seqs )
            if( buf[7] == '+' ):
                s_pos = int( buf[9] )
                e_pos = int( buf[10] )
                j = common.detect_index( buf[6], s_pos, e_pos, seqs )
            else:
                s_pos = int( buf[8] ) - int( buf[9] )
                e_pos = int( buf[8] ) - int( buf[10] )
                j = common.detect_index( buf[6], e_pos, s_pos, seqs )
            color = heatmap.convert_identity2color( float( buf[12][:-1] ))

            if( i == -1 or j == -1 ):
                continue
            
            x, y = convert_position2coord( seqs[i][buf[1]], seqs[j][buf[6]], int( buf[4] ), int( buf[5] ), s_pos, e_pos, size.margin_bw_scaffold_alignment, size.histograms[i], size.histograms[j] )
            ax.fill( x, y, color=color, alpha=heatmap.alpha, lw=0 )


def plot_alignment4mummer( seqs, ax, heatmap, size, fn ):
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( )
            i = common.detect_index( buf[11], int( buf[0] ), int( buf[1] ), seqs )
            j = common.detect_index( buf[12], int( buf[3] ), int( buf[4] ), seqs )

            color = heatmap.convert_identity2color( float( buf[9] ))

            if( i == -1 or j == -1 ):
                continue

            x, y = convert_position2coord( seqs[i][buf[11]], seqs[j][buf[12]], int( buf[0] ), int( buf[1] ), int( buf[3] ), int( buf[4] ), size.margin_bw_scaffold_alignment, size.histograms[i], size.histograms[j] )
            ax.fill( x, y, color=color, alpha=heatmap.alpha, lw=0 )

            
def cal_min_identity4original( fn ):
    min_identity = 0
    BIN = 10
    dict = {}
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( '\t' )
            dict[ int( float( buf[6] )/BIN ) * BIN ] = 0

    min_identity = min( dict )
    return min_identity


def cal_min_identity4blastn( fn ):
    min_identity = 0
    BIN = 10
    dict = {}
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( '\t' )
            dict[ int( float( buf[2] )/BIN ) * BIN ] = 0

    min_identity = min( dict )
    return min_identity


def cal_min_identity4lastz( fn ):
    min_identity = 0
    BIN = 10
    dict = {}
    with open( fn , 'r' ) as file:
        for line in file:
            if( line[0:1] == "#" ):
                continue
            buf = line.rstrip( '\n' ).split( '\t' )
            dict[ int( float( buf[12][:-1] )/BIN ) * BIN ] = 0

    min_identity = min( dict )
    return min_identity


def cal_min_identity4mummer( fn ):
    min_identity = 0
    BIN = 10
    dict = {}
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( )
            dict[ int( float( buf[9] )/BIN ) * BIN ] = 0

    min_identity = min( dict )
    return min_identity
