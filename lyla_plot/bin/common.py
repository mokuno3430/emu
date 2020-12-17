import argparse
from matplotlib import colors as mcolors
import os.path

def get_args():
    parser = argparse.ArgumentParser( formatter_class=argparse.MetavarTypeHelpFormatter )
    
    parser.add_argument( '-i', '--input', help='scaffold_info.tsv (arranged from bottom to top)', type=str, nargs='*', required=True )
    parser.add_argument( '-a', '--alignment', help='alignment.tsv', type=str, nargs='*' )
    parser.add_argument( '--blastn', help='blastn.tsv (-outfmt 6)', type=str, nargs='*' )
    parser.add_argument( '--lastz', help='lastz.tsv (--format=general)', type=str, nargs='*' )
    parser.add_argument( '--mummer', help='show-coords.tsv (--format=show-coords -H)', type=str, nargs='*' )
    parser.add_argument( '--gff3', help='gff3 file', type=str, nargs='*' )
    parser.add_argument( '--hist', help='histogram file (scaffold_ID\tposition\tvalue)', type=str, nargs='*' )
    parser.add_argument( '--mark_v', help='fileformat in the air', type=str )
    parser.add_argument( '--out', help='optional: prefix of pdf file (default out)', type=str, default='out' )
    parser.add_argument('--scaffold_layout', help='optional: (default center)', choices=['left', 'center', 'right'], default='center', type=str )
    parser.add_argument( '--min_identity', help='optional: (default -1, -1 means auto)', type=int, default=-1 )
    parser.add_argument( '--max_identity', help='optional: (default 100)', type=int, default=100 )
    parser.add_argument( '--margin_bw_scaffolds', help='optional: (default -1, -1 means auto)', type=float, default=-1 )
    parser.add_argument( '--alignment_height', help='optional: (default 1.5)', type=float, default=1.5 )
    parser.add_argument( '--xlim_max', help='optional: (default -1, -1 means auto)', type=float, default=-1 )
    parser.add_argument( '--scaffold_font_size', help='optional: (default 0, 0 means not shown)', type=float, default=0 )
    parser.add_argument('--colormap', help='optional: colormap for identity of alignments ( 0:binary, 1:bone_r, 2:inferno_r, 3:hot_r, 4:YlGnBu, 5:original) (default 0)', choices=[ 0, 1, 2, 3, 4, 5 ], default=0, type=int )
    args = parser.parse_args()

    return ( args )


class Size:
    def __init__( self, seqs, margin_bw_scaffolds, xlim_max, alignment_height ):
        length = cal_total_length( seqs )        
        self.alignment_height = alignment_height
        self.alignments = [ alignment_height ] * len( seqs )
        self.alignments[-1] = 0
        self.histogram_height = 1.0
        self.histograms = [ 0 ] * len( seqs )
        self.margin_bw_scaffolds = int( max( length )/200 ) if margin_bw_scaffolds == -1 else margin_bw_scaffolds
        self.xlim_max = self.cal_xlim_max( length, seqs, xlim_max )
        self.bottom_margin = 1.0
        self.top_margin = 0.5
        self.ylim_max = 5
        self.figsize_inch = [8,4]
        self.margin_bw_scaffold_alignment = 0.03
        
    def cal_xlim_max( self, length, seqs, xlim_max ):
        sum_len = []
        for i in range( len( length ) ):
            sum_len.append( length[i] + self.margin_bw_scaffolds * ( len( seqs[i] ) - 1 ) )
        return max( sum_len ) if max( sum_len ) > xlim_max else xlim_max
    
    def set_scaffold_layout( self, seqs, layout ):
        length = cal_total_length( seqs )
        for i in range( len( length ) ):
            length[i] += self.margin_bw_scaffolds * ( len( seqs[i] ) - 1 )

        y = self.bottom_margin
        for i in range( len( seqs )):
            x = self.__set_xvalue( layout, length[i] )            
            for scaf in seqs[i]:
                seqs[i][scaf].set_origins( x, y )
                x += seqs[i][scaf].length + self.margin_bw_scaffolds
            y += self.alignments[i] + self.histograms[i]
        y += self.top_margin
        self.ylim_max = y
        
    def __set_xvalue( self, layout, length ):
        if( layout == 'center' ):
            return int( ( self.xlim_max - length )/2 )
        if( layout == 'left' ):
            return 0
        return self.xlim_max - length
        
    def set_histogram_space( self, seqs, files ):
        if files is None: 
            return
        for fn in files:
            if not os.path.isfile( fn ):
                continue
            with open( fn , 'r' ) as file:
                for line in file:
                    buf = line.rstrip( '\n' ).split( '\t' )
                    i = detect_index( buf[0], int( buf[2] ), int( buf[2] ), seqs )
                    if( i == -1 ):
                        continue
                    self.histograms[ i ] = self.histogram_height
    

    def output_parameters( self ):
        print( '\n##Size paramenters:' )
        print( '  figsize_inch: %.2f,%.2f' % ( self.figsize_inch[0], self.figsize_inch[1] ))
        print( '  xlim_max: %d' % ( self.xlim_max ))
        print( '  ylim_max: %.2f' % ( self.ylim_max ))
        print( '  margin_bw_scaffolds: %d' % ( self.margin_bw_scaffolds ))
        print( '  alignment_height: %.2f' % ( self.alignment_height ))
        print( '  histogram_height: %.2f' % ( self.histogram_height ))
        print( '  top_margin: %.2f' % ( self.top_margin ))
        print( '  bottom_margin: %.2f' % ( self.bottom_margin ))
        print( '' )


class Text:
    def __init__( self, label, size, x, y, ha, va, color='black' ):
        self.label = label
        self.size = size
        self.origin_x = x
        self.origin_y = y
        self.color = color
        self.ha = ha
        self.va = va        

    def output( self, ax ):
        if self.label == 'BLANK':
            return
        if self.size == 0:
            return
        ax.text( self.origin_x, self.origin_y, self.label, fontsize = self.size, color = self.color, ha=self.ha, va=self.va )

        
class Polygon:
    def __init__( self, x, y, color, alpha ):
        self.x = x
        self.y = y
        self.color = Color( color, alpha )
        
    def plot( self, ax ):
        ax.fill( self.x, self.y, color=self.color.color, alpha=self.color.alpha, lw=0 )
        
        
class Color:
    replaced_color = 'grey'
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

    def __init__( self, color, alpha ):
        self.color = self.__evaluate_input( color )
        self.alpha = alpha
        
    def __evaluate_input( self, color ):
        if color in Color.colors:
            return color
        if color[0:1] == '#' and len( color ) == 7:
            RGB = (int(color[1:3],16),int(color[3:5],16),int(color[5:7],16))
            if( RGB[0] <= 255 and RGB[1] <= 255 and RGB[2] <= 255 ):
                return color
        print( "invalid color :{c}".format( c = color ) )
        return Color.replaced_color

        
def cal_total_length( seqs ):
    length = []
    for sample in seqs:
        total = 0
        for scaf in sample:
            total += sample[scaf].length
        length.append( total )
    return length


def detect_index( ID, start, end, seqs ):
    for j in range( len( seqs )):
        if ID in seqs[j].keys():
            if( seqs[j][ID].name == 'BLANK' ):
                return -1
            if( seqs[j][ID].start <= start and end <= seqs[j][ID].end ):
                return j
    return -1
