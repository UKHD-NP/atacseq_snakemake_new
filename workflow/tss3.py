import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pybedtools
import metaseq
import logging
import os
from matplotlib import pyplot as plt

def make_tss_plot(
    bam_file,
    tss,
    chromsizes,
    read_len,
    out_plot,
    out_plot_large,
    sample_name=None,
    bins=400,
    bp_edge=2000,
    processes=8,
    greenleaf_norm=True
):
    logging.info('Generating tss plot...')
    
    # Extract sample name if not provided
    if sample_name is None:
        sample_name = '_'.join(os.path.basename(bam_file).split('_')[:3])  # Extract from BAM filename
        if not sample_name:
            sample_name = "Sample"
    
    # Load the TSS file
    tss_obj = pybedtools.BedTool(tss)
    tss_ext = tss_obj.slop(b=bp_edge, g=chromsizes)
    
    # Validate extended regions (basic check)
    if len(tss_ext) == 0:
        logging.error("No valid TSS regions after extension")
        sys.exit(1)
    
    # Load the BAM file
    bam = metaseq.genomic_signal(bam_file, 'bam')
    
    # Use standard ATAC-seq offset 
    shift_width = -read_len // 2
    
    bam_array = bam.array(
        tss_ext,
        bins=bins,
        shift_width=shift_width,
        processes=processes,
        stranded=True
    )
    
    # Initialize enrichment score
    tss_enrichment = 1.0
    
    if greenleaf_norm:
        # Ensure minimum number of edge bins (at least 10% of total bins or 20, whichever is smaller)
        min_edge_bins = min(20, bins // 10)
        calculated_edge_bins = int(100 / (2 * bp_edge / bins))
        num_edge_bins = max(min_edge_bins, calculated_edge_bins)
        
        # Ensure edge bins don't exceed half the total bins
        num_edge_bins = min(num_edge_bins, bins // 4)
        
        bin_means = bam_array.mean(axis=0)
        tss_enrichment = calculate_tss_enrichment(bin_means, bp_edge, bins)
        logging.info("TSS enrichment score: {:.2f}".format(tss_enrichment))
        
        with open(out_plot.replace('.pdf','.score.txt'), 'w') as f:
            f.write(str(tss_enrichment) + '\n')
        
        # Calculate average noise with safeguards
        edge_left = bin_means[:num_edge_bins]
        edge_right = bin_means[-num_edge_bins:]
        
        # Filter out extreme outliers that might skew normalization
        edge_values = np.concatenate([edge_left, edge_right])
        edge_values = edge_values[edge_values > 0]  # Remove zero values
        
        if len(edge_values) > 0:
            avg_noise = np.mean(edge_values)
            
            # Apply minimum threshold to prevent division by very small numbers
            min_noise_threshold = np.percentile(bin_means[bin_means > 0], 5) if len(bin_means[bin_means > 0]) > 0 else 1e-6
            avg_noise = max(avg_noise, min_noise_threshold)
        else:
            logging.warning("No valid edge values found, using fallback normalization")
            avg_noise = np.median(bin_means[bin_means > 0]) if len(bin_means[bin_means > 0]) > 0 else 1.0
        
        logging.info("Average noise level: {:.6f}".format(avg_noise))
        
        # Apply normalization 
        normalized_array = bam_array / avg_noise
        
        bam_array = normalized_array
    else:
        # Standard RPM normalization
        mapped_reads = bam.mapped_read_count()
        if mapped_reads == 0:
            logging.error("No mapped reads found in BAM file")
            sys.exit(1)
        bam_array = bam_array / (mapped_reads / 1e6)
    
    # Generate a line plot
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    x = np.linspace(-bp_edge, bp_edge, bins)
    
    mean_profile = bam_array.mean(axis=0)
    ax.plot(x, mean_profile, color='r', label='Mean', linewidth=1.5)
    ax.axvline(0, linestyle=':', color='k', alpha=0.7)
    
    ax.set_xlabel('Distance from TSS (bp)')
    if greenleaf_norm:
        ax.set_ylabel('TSS enrichment score')
    else:
        ax.set_ylabel('Average read coverage (per million mapped reads)')
    
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Add sample name and enrichment score to the plot
    y_pos_high = ax.get_ylim()[1] * 0.95  # Position near top of plot
    y_pos_mid = ax.get_ylim()[1] * 0.88   # Slightly lower position
    
    # Add sample name (top-left)
    ax.text(0.02, 0.98, 'Sample: {}'.format(sample_name), 
            transform=ax.transAxes, 
            fontsize=12, 
            fontweight='bold',
            verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    # Add enrichment score (top-right)
    if greenleaf_norm:
        ax.text(0.98, 0.88, 'TSS Enrichment: {:.2f}'.format(tss_enrichment), 
                transform=ax.transAxes, 
                fontsize=12, 
                fontweight='bold',
                horizontalalignment='right',
                verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.8))
    else:
        # For RPM normalization, show total mapped reads instead
        mapped_reads = bam.mapped_read_count()
        ax.text(0.88, 0.88, 'Mapped reads: {:.1f}M'.format(mapped_reads/1e6), 
                transform=ax.transAxes, 
                fontsize=12, 
                fontweight='bold',
                horizontalalignment='right',
                verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.8))
    
    # Add additional info (bottom-right)
    info_text = 'Bins: {} | Window: +-{}bp'.format(bins, bp_edge)
    ax.text(0.98, 0.02, info_text, 
            transform=ax.transAxes, 
            fontsize=10, 
            horizontalalignment='right',
            verticalalignment='bottom',
            style='italic',
            alpha=0.7)
    
    fig.savefig(out_plot, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # Heatmap-like plot (also add sample name here)
    upper_prct = 99
    array_values = bam_array.ravel()
    array_values = array_values[array_values > 0]  # Remove zeros for percentile calculation
    
    if len(array_values) == 0 or np.percentile(array_values, upper_prct) == 0.0:
        upper_prct = 100.0
    
    plt.rcParams['font.size'] = 8
    fig = metaseq.plotutils.imshow(
        bam_array,
        x=x,
        figsize=(5, 10),
        vmin=5, vmax=upper_prct, percentile=True,
        line_kwargs=dict(color='k', label='All'),
        fill_kwargs=dict(color='k', alpha=0.3),
        sort_by=bam_array.mean(axis=1)
    )
    
    # Add sample name to heatmap plot
    fig.suptitle('{} - TSS Enrichment Heatmap'.format(sample_name), fontsize=14, fontweight='bold', y=0.98)
    if greenleaf_norm:
        # Add enrichment score as subtitle
        fig.text(0.5, 0.94, 'TSS Enrichment Score: {:.2f}'.format(tss_enrichment), 
                ha='center', va='top', fontsize=12,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.8))
    
    fig.savefig(out_plot_large, dpi=300, bbox_inches='tight')
    plt.close(fig)

def calculate_tss_enrichment(profile, bp_edge, bins):
    """Calculate TSS enrichment score with improved robustness"""
    bin_size = (2 * bp_edge) / bins
    center_window = max(1, int(100 / bin_size))  # +-100bp, ensure at least 1 bin
    edge_window = max(1, int(100 / bin_size))    # +-100bp, ensure at least 1 bin
    
    # Ensure windows don't exceed available bins
    center_window = min(center_window, bins // 4)
    edge_window = min(edge_window, bins // 4)
    
    center_start = bins // 2 - center_window // 2
    center_end = bins // 2 + center_window // 2
    
    # Ensure valid indices
    center_start = max(0, center_start)
    center_end = min(bins, center_end)
    edge_window = min(edge_window, center_start)  # Don't overlap with center
    
    # Edge regions: first and last bins
    if edge_window > 0:
        center_signal = profile[center_start:center_end].mean()
        edge_signal = np.concatenate([profile[:edge_window], profile[-edge_window:]]).mean()
        
        # Avoid division by zero
        if edge_signal > 0:
            enrichment = center_signal / edge_signal
        else:
            enrichment = 1.0  # Default when no edge signal
    else:
        enrichment = 1.0
    
    return max(enrichment, 0.0)  # Ensure non-negative enrichment


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Generate TSS plots.")
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--tss', required=True, help='TSS BED file')
    parser.add_argument('--chromsizes', required=True, help='Chromosome sizes file')
    parser.add_argument('--read_len', type=int, required=True, help='Read length')
    parser.add_argument('--out_plot', required=True, help='Output plot file')
    parser.add_argument('--out_plot_large', required=True, help='Output large plot file')
    parser.add_argument('--sample_name', help='Sample name to display on plot')
    parser.add_argument('--bins', type=int, default=400, help='Number of bins')
    parser.add_argument('--bp_edge', type=int, default=2000, help='Base pairs from TSS edge')
    parser.add_argument('--processes', type=int, default=8, help='Number of processes')
    parser.add_argument('--greenleaf_norm', action='store_true', help='Use Greenleaf normalization')
    
    args = parser.parse_args()
    
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s"
    )
    
    try:
        make_tss_plot(
            bam_file=args.bam,
            tss=args.tss,
            chromsizes=args.chromsizes,
            read_len=args.read_len,
            out_plot=args.out_plot,
            out_plot_large=args.out_plot_large,
            sample_name=args.sample_name,
            bins=args.bins,
            bp_edge=args.bp_edge,
            processes=args.processes,
            greenleaf_norm=args.greenleaf_norm
        )
        logging.info("TSS plot generation completed successfully")
    except Exception as e:
        logging.exception("Exception occurred: {}".format(e))
        sys.exit(1)