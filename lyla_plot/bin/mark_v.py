import matplotlib.pyplot as plt
import sys
import csv
import common

class Vmark:

    def __init__( self, scaffold, buf, width, height, color ):
        x1, y1 = self.__convert_position2coord( scaffold, buf ) 
        self.width = width
        self.height = height
        self.x = [ x1, x1 + self.width /2, x1 - self.width /2]
        self.y = [ y1, y1 + self.height, y1 + self.height]
        #self.y = [ y1 + 1, y1 + self.height + 1, y1 + self.height + 1]
        self.color = common.Color( color, 1 )
        
    def __convert_position2coord( self, scaffold, buf ):
        x1 = scaffold.convert_position2xcoord( ( int( buf[3] ) + int( buf[4] ) ) / 2 )
        y1 = scaffold.convert_position2ycoord( scaffold.height )
        return x1, y1

    def plot( self, ax ):
        ax.fill( self.x, self.y, color=self.color.color, alpha=self.color.alpha, lw=0, zorder=3 )
    

def plot_mark_v( seqs, ax, size, fn ):
    height = size.ylim_max /50
    width = size.xlim_max / 150
    #color = 'tomato'
    color = '#043c78'
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( '\t' )
            i = common.detect_index( buf[0], int( buf[3] ), int( buf[4] ), seqs )

            if( i == -1 ):
                continue
            aninstance = Vmark( seqs[i][buf[0]], buf, width, height, color )
            aninstance.plot( ax )
            
