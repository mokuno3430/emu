import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os.path

import common
import scaffold as scf
import alignment
import mark_v
import gff
import histogram

argv = sys.argv

def func_set_axes( ax, size ):
    ##set data range
    ax.set_xlim(0, size.xlim_max)
    ax.set_ylim(0, size.ylim_max)
                
    ##non-display axis
    ax.spines["right"].set_color("none")  
    ax.spines["left"].set_color("none")   
    ax.spines["top"].set_color("none")    
    ax.spines["bottom"].set_color("none") 
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.tick_params( length=0 )

    ax.patch.set_alpha(0.0)

    
def func_print_messages():
    print( 'start' ) 
    print( ' '.join( argv ))

    
def main():

    args = common.get_args()
    func_print_messages( )
    seqs = scf.input_scaffoldtsv( args.input )
    size=common.Size( seqs, args.margin_bw_scaffolds, args.xlim_max, args.alignment_height )
    #histograms = histogram.set_space( args.hist, seqs, size.histogram_height )
    size.set_histogram_space( seqs, args.hist ) 
    size.set_scaffold_layout( seqs, args.scaffold_layout )
    size.output_parameters()
        
    fig = plt.figure( figsize=size.figsize_inch )
    ax = fig.add_subplot(111)
    fig.patch.set_alpha( 0.0 )
    func_set_axes( ax, size )
    
    scf.plot_scaffolds( ax, seqs, args.scaffold_font_size )
    
    ##plot scale bar
    scalebar = scf.Scalebar( size )
    scalebar.plot( ax )
    scalebar.output_parameters()
    
    ##plot alignment
    max_identity = args.max_identity
    input_formats = [ args.alignment, args.blastn, args.lastz, args.mummer ]
    func_plot_alignmment = [ alignment.plot_alignment4original, alignment.plot_alignment4blastn, alignment.plot_alignment4lastz, alignment.plot_alignment4mummer ] 
    valid_files = alignment.count_alignment_files( args )
    if valid_files == 0:
        pass
    else:
        min_identity = alignment.set_min_identity( args )
        ##set colormap
        heatmap = alignment.Colormap( min_identity, max_identity, args.colormap )
        heatmap.output_parameters()
        ##set and plot colormap legend
        heatmap_legend = alignment.Colorbox( size )
        heatmap_legend.plot( ax, heatmap )
        heatmap_legend.output_parameters()

        for files, func_plot in zip( input_formats, func_plot_alignmment ):
            if files is None:
                continue
            for fn in files:
                if not os.path.isfile( fn ):
                    continue
                func_plot( seqs, ax, heatmap, size, fn )
                     
    ##plot mark_v
    if args.mark_v is not None: 
        if os.path.isfile( args.mark_v ):
            mark_v.plot_mark_v( seqs, ax, size, args.mark_v )

    ##plot gene
    if args.gff3 is not None: 
        for fn in args.gff3:
            if not os.path.isfile( fn ):
                continue
            gff.plot_genes( seqs, ax, size, fn )

    ##plot histogram
    histogram.plot_background( seqs, ax, size )
    if args.hist is not None: 
        for fn in args.hist:
            if not os.path.isfile( fn ):
                continue
            histogram.plot_histogram( seqs, ax, size, fn )

    pdf_file = args.out + '.pdf'
    pp = PdfPages( pdf_file )
    pp.savefig( fig, bbox_inches='tight' )
    pp.close()
    

if __name__ == '__main__':
    main()
