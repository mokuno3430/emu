import matplotlib.pyplot as plt
import sys
import csv
import common

class Gene:

    def __init__( self, scaffold, buf, height, tail_length, color ):
        self.start, self.end, self.strand = self.__convert_position2coord( scaffold, buf )
        self.y_mid_coord = scaffold.convert_position2ycoord( scaffold.height /2 )
        self.width = self.end - self.start
        self.height = height
        self.color = common.Color( color, 1 )
        self.x, self.y  = self.__set_gene_coord( tail_length )

    def __set_gene_coord( self, tail_length ):
        x1 = self.start
        x3 = self.end
        y1 = self.y_mid_coord - self.height / 2
        #y1 = self.y_mid_coord - self.height * 0.8
        y2 = self.y_mid_coord
        y3 = self.y_mid_coord + self.height / 2
        #y3 = self.y_mid_coord + self.height * 0.8
        if( self.strand == '+' ): # => or >
            x2 = self.end - tail_length if self.width > tail_length else self.start            
            x = [ x1, x2, x3, x2, x1 ]
            y = [ y1, y1, y2, y3, y3 ]
        elif( self.strand == '-' ): # <= or <
            x2 = self.start + tail_length if self.width > tail_length else self.end
            x = [ x1, x2, x3, x3, x2 ]
            y = [ y2, y1, y1, y3, y3 ]
        else: 
            x = [ x1, x3, x3, x1 ]
            y = [ y1, y1, y3, y3 ]
        return x, y

    def __convert_position2coord( self, scaffold, buf ):
        if( scaffold.strand == '+' ):
            start = scaffold.convert_position2xcoord( int( buf[3] ) )
            end = scaffold.convert_position2xcoord( int( buf[4] ) )
            return start, end, buf[6]
        else:
            start = scaffold.convert_position2xcoord( int( buf[4] ) )
            end = scaffold.convert_position2xcoord( int( buf[3] ) )
            if( buf[6] == '+' ): #逆にする
                return  start, end, '-'
            return start, end, '+'

    def plot( self, ax ):
        ax.fill( self.x, self.y, color=self.color.color, alpha=self.color.alpha, lw=0 )
        
    
def plot_genes( seqs, ax, size, fn ):
    color = 'darkorange'
    RATIO = 3
    #color = 'tomato'
    
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( '\t' )
            if( len( buf ) < 4 ):
                continue
            if( buf[2] != 'gene' ):
                continue
            
            i = common.detect_index( buf[0], int( buf[3] ), int( buf[4] ), seqs )
            if( i == -1 ):
                continue
            
            height = seqs[i][buf[0]].height * RATIO
            tail_length = ( height / 2 ) * ( size.xlim_max / size.ylim_max ) * ( size.figsize_inch[1] / size.figsize_inch[0] )
            aninstance = Gene( seqs[i][buf[0]], buf, height, tail_length, color )
            aninstance.plot( ax )
            
