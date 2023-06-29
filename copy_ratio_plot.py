"""
param filename, a .tsv file containing functional annotations
param chrom, the chromosome to plot
param start, the genomic coordinate at which to start plotting
param end, the genomic coordinate at which to finish plotting
return, a dictionary mapping plotted intervals to GC content
return, a dictionary mapping plotted intervals to segdup content
return, a dictionary mapping plotted intervals to mappability 
"""
def read_functional_annotations(filename, chrom, start, end):
    gc_dict = dict()
    sd_dict = dict()
    m_dict = dict()
    with open(filename, "r") as f:
        header = f.readline().split()
        contig_idx = None
        istart_idx = None
        iend_idx = None
        gc_idx = None
        sd_idx = None
        m_idx = None
        for i in range(len(header)):
            if header[i] == "CONTIG":
                contig_idx = i
            elif header[i] == "START":
                istart_idx = i
            elif header[i] == "END":
                iend_idx = i
            elif header[i] == "GC_CONTENT":
                gc_idx = i
            elif header[i] == "MAPPABILITY":
                m_idx = i
            elif header[i] == "SEGMENTAL_DUPLICATION_CONTENT":
                sd_idx = i
        line = f.readline()
        while line != "":
            spline = line.split()
            if spline[contig_idx] == chrom and int(spline[iend_idx]) >= start and int(spline[istart_idx]) <= end:
                interval = (int(spline[istart_idx]), int(spline[iend_idx]))
                gc_dict[interval] = float(spline[gc_idx])
                sd_dict[interval] = float(spline[sd_idx])
                m_dict[interval] = float(spline[m_idx])
            line = f.readline()
    return gc_dict, sd_dict, m_dict

"""
param filename, a tabix-indexed .tsv file containing copy ratio data
param chrom, the chromosome to plot
param start, the genomic coordinate at which to start plotting
param end, the genomic coordinate at which to finish plotting
return, an intervals x samples array containing copy ratio data
return, an intervals x 3 array specifying the boundaries of each interval as (chrom, start, end)
return, an array of sample names
"""
def read_copy_ratios(filename, chrom, start, end):
    f = pysam.TabixFile(filename)
    h = f.header[0].strip("#").split()
    contig_idx = None
    istart_idx = None
    iend_idx = None
    for i in range(len(h)):
        if h[i] == "CONTIG":
            contig_idx = i
        elif h[i] == "START":
            istart_idx = i
        elif h[i] == "END":
            iend_idx = i
    if end != float("inf"):
        row_iter = f.fetch(reference=chrom, start=start, end=end)
    else:
        row_iter = f.fetch(reference=chrom, start=start)
    cr_list = list()
    cr_region_list = list()
    sample_names = np.delete(h, [contig_idx, istart_idx, iend_idx])
    for row in row_iter:
        r = row.split()
        cr_list.append(np.delete(r, [contig_idx, istart_idx, iend_idx]))
        cr_region_list.append([r[contig_idx], int(r[istart_idx]), int(r[iend_idx])])
    return np.array(cr_list, dtype="float64"), np.array(cr_region_list, dtype="object"), sample_names

"""
Ensures that displayed sample names do not overlap each other
param sample_labels, a dictionary mapping y-coordinates to sample names
parram y_buffer, the amount of vertical space required between sample names
param max_iters, the maximum number of successive iterations for clustering
return, a dictionary mapping cluster y-coordinate intervals to tuples of original y-coordinates
"""
def cluster_labels(sample_labels, y_buffer, max_iters=10):
    label_clusters = dict() # (min_y, max_y) after buffer -> tuple of original y's
    for l in sample_labels:
        label_clusters[(l - y_buffer / 2, l + y_buffer / 2)] = (l,)
    changed = True
    i = 0
    while changed and i < max_iters:
        changed = False
        i += 1
        to_merge = dict() # original cluster -> new cluster
        clusters = list(label_clusters.keys())
        for j1 in range(len(clusters)):
            for j2 in range(j1 + 1, len(clusters)):
                l1, l2 = (clusters[j1], clusters[j2])
                if l1[0] < l2[1] and l2[0] < l1[1]: # overlap
                    new_cluster = set((l1, l2))
                    if l1 in to_merge:
                        new_cluster.update(to_merge[l1])
                    if l2 in to_merge:
                        new_cluster.update(to_merge[l2])
                    fully_updated = False
                    while not fully_updated:
                        fully_updated = True
                        for l in new_cluster:
                            if l in to_merge:
                                for m in to_merge[l]:
                                    if m not in new_cluster:
                                        new_cluster.add(m)
                                        fully_updated = False
                    new_cluster = tuple(new_cluster)
                    for l in new_cluster:
                        to_merge[l] = new_cluster
        if len(to_merge) > 0:
            changed = True
            added = set()
            for l in to_merge:
                new_cluster = to_merge[l]
                if l not in added:
                    curr_min = float("inf")
                    curr_max = float("-inf")
                    labels = 0
                    contents = list()
                    for m_min, m_max in new_cluster:
                        curr_min = min(curr_min, m_min)
                        curr_max = max(curr_max, m_max)
                        labels += (m_max - m_min) / y_buffer
                        contents.extend(label_clusters[(m_min, m_max)])
                        if (m_min, m_max) in label_clusters:
                            label_clusters.pop((m_min, m_max))
                            added.add((m_min, m_max))
                    size_incr = labels * y_buffer - (curr_max - curr_min)
                    new_max = curr_max + size_incr / 2
                    new_min = curr_min - size_incr / 2
                    label_clusters[(new_min, new_max)] = tuple(contents)
    return label_clusters

"""
param args, a namespace object containing input parameters
"""
def create_plot(args):
    
    # Read in data
    cr_matrix, cr_regions, sample_names = read_copy_ratios(args.data, args.chrom, args.start, args.end)
    if args.functional_annotations is not None:
        gc_dict, sd_dict, m_dict = read_functional_annotations(args.functional_annotations, args.chrom, args.start, args.end)
    
    # Basic parameters and plot panes
    n_regions = cr_matrix.shape[0]
    n_samples = cr_matrix.shape[1]
    to_show = [args.show_missing, args.show_gc, args.show_segdup, args.show_map]
    if sum(to_show) > 0:
        subplot = True
        subplot_ratio = 10 / (sum(to_show))
        subplot_labels = np.array(["Missing fraction", "GC content", "Segdup content", "Mappability"])
        subplot_labels = subplot_labels[np.where(to_show)]
    else:
        subplot = False
    
    # Generate bin data
    n_bins = min(args.max_bins, n_regions)
    if args.max_bins >= n_regions:
        bin_matrix = cr_matrix
        bin_regions = cr_regions
        bin_idx_regions = np.array([[i,i+1] for i in range(n_regions)])
    else:
        bin_matrix = np.zeros((args.max_bins, n_samples))
        bin_regions = np.zeros((args.max_bins, 3), dtype="object")
        bin_idx_regions = np.zeros((args.max_bins, 2))
        if args.evenly_pooled_bins:
            bin_size = n_regions // args.max_bins
            bin_size_remainder = n_regions % args.max_bins
            i_bin_start = 0
            for i in range(args.max_bins):
                if i < bin_size_remainder:
                    i_bin_size = bin_size + 1
                else:
                    i_bin_size = bin_size
                i_bin_end = i_bin_start + i_bin_size - 1
                bin_regions[i,:] = [cr_regions[i_bin_start,0], cr_regions[i_bin_start,1], cr_regions[i_bin_end,2]]
                if args.bin_median:
                    bin_matrix[i,:] = np.nanmedian(cr_matrix[i_bin_start:i_bin_end+1,:], axis=0)
                else:
                    bin_matrix[i,:] = np.nanmean(cr_matrix[i_bin_start:i_bin_end+1,:], axis=0)
                bin_idx_regions[i,:] = [i_bin_start, i_bin_end+1]
                i_bin_start = i_bin_end + 1
        else:
            bin_size = (cr_regions[-1,2] - cr_regions[0,1]) // args.max_bins
            bin_size_remainder = (cr_regions[-1,2] - cr_regions[0,1]) % args.max_bins
            i_bin_start = cr_regions[0,1]
            r = 0
            for i in range(args.max_bins):
                if i < bin_size_remainder:
                    i_bin_size = bin_size + 1
                else:
                    i_bin_size = bin_size
                i_bin_end = i_bin_start + i_bin_size
                bin_regions[i,:] = [cr_regions[r,0], i_bin_start, i_bin_end]
                bin_regs = list()
                r_center = np.mean(cr_regions[r,1:])
                while r < n_regions and (r_center >= i_bin_start and r_center < i_bin_end):
                    bin_regs.append(r)
                    r += 1
                    if r < n_regions:
                        r_center = np.mean(cr_regions[r,1:])
                if args.bin_median:
                    bin_matrix[i,:] = np.nanmedian(cr_matrix[bin_regs,:], axis=0)
                else:
                    bin_matrix[i,:] = np.nanmean(cr_matrix[bin_regs,:], axis=0)
                if len(bin_regs) > 0:
                    bin_idx_regions[i,:] = [min(bin_regs), max(bin_regs)+1]
                else:
                    bin_idx_regions[i,:] = [float("nan"), float("nan")]
                i_bin_start = i_bin_end
            
    # Get points to plot
    plot_points = list()
    to_skip = set()
    for r in range(n_bins):
        if np.all(np.isnan(bin_matrix[r,:])):
            to_skip.add(r)
    for s in range(n_samples):
        sample_list = list()
        for r in range(n_bins):
            if r not in to_skip:
                if args.evenly_spaced_intervals:
                    sample_list.append((bin_idx_regions[r,0], bin_matrix[r,s]))
                    sample_list.append((bin_idx_regions[r,1], bin_matrix[r,s]))
                else:
                    sample_list.append((bin_regions[r,1], bin_matrix[r,s]))
                    sample_list.append((bin_regions[r,2], bin_matrix[r,s]))
        plot_points.append(sample_list)
    
    # Create plot structure
    if subplot:
        fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [subplot_ratio, 1]}, figsize=(6.4, 4.8))
        fig.subplots_adjust(hspace=0)
        ax2.set_yticks([])
    else:
        fig = plt.figure(figsize=(6.4, 4.8))
        ax = plt.gca()
    if args.evenly_spaced_intervals:
        (x_min, x_max) = (0, n_regions)
    else:
        (x_min, x_max) = (np.min(cr_regions[:,1:]), np.max(cr_regions[:,1:]))
    x_len = x_max - x_min
    x_buffer = x_len * 0.05
    (x_min, x_max) = (x_min - x_buffer, x_max + x_buffer)
    (y_min, y_max) = (0, 5.25)
    
    # Display interval locations and plot copy ratio data
    if args.show_intervals:
        for r in range(n_regions):
            if args.evenly_spaced_intervals:
                ax.add_patch(Rectangle((r, y_min), 1 , y_max - y_min, color="0.9", zorder=0))
            else:
              ax.add_patch(Rectangle((cr_regions[r,1], y_min), cr_regions[r,2] - cr_regions[r,1] , y_max - y_min, color="0.9", zorder=0))
        sample_color = "0.65"
    else:
        sample_color = "0.8"
    if args.show_samples:
        sample_labels = dict() # y coordinate -> sample names
        label_line_bin = dict() # y coordinate -> x-coordinate of last non-missing point for sample
    for s in range(n_samples):
        label = False
        if args.auto_highlight:
            sample_dev_sign = np.sign(bin_matrix[:,s] - args.exp_cr)
            sample_aberr = np.abs(bin_matrix[:,s] - args.exp_cr) > args.aberrant_threshold
            sample_adj_aberr = [True if sample_aberr[i] and sample_aberr[i+1] else False for i in range(n_bins-1)]
            sample_adj_aberr = [sample_adj_aberr[i] if sample_dev_sign[i] == sample_dev_sign[i+1] else False for i in range(n_bins-1)]
        if sample_names[s] in args.sample_highlight:
            color = args.sample_color
            zorder = 3
            label = True
        elif args.auto_highlight and np.any(sample_adj_aberr):
            aber_gain_frac = np.mean(bin_matrix[:,s][np.where(np.abs(bin_matrix[:,s] - args.exp_cr) > args.aberrant_threshold)] > args.exp_cr)
            sample_mean = np.mean(bin_matrix[:,s])
            if aber_gain_frac > 0.5 or (aber_gain_frac == 0.5 and sample_mean > 0.5):
                color = args.aberrant_dup_color
            elif aber_gain_frac < 0.5 or (aber_gain_frac == 0.5 and sample_mean < 0.5):
                color = args.aberrant_del_color
            else:
                color = np.random.choice([args.aberrant_dup_color, args.aberrant_del_color])
            zorder = 2
            label = True
        else:
            color = sample_color
            zorder = 1
        if args.show_samples and label:
            i = len(bin_matrix) - 1
            y = bin_matrix[i,s]
            while i >= 0 and np.isnan(y):
                i -= 1
                y = bin_matrix[i,s]
            if y not in sample_labels:
                sample_labels[y] = sample_names[s]
                if args.evenly_spaced_intervals:
                    label_line_bin[y] = (bin_idx_regions[i,1],)
                else:
                    label_line_bin[y] = (bin_regions[i,2],)
            else:
                sample_labels[y] = (sample_labels[y], sample_names[s])
                if args.evely_spaced_intervals:
                    label_line_bin[y] = tuple(list(label_line_bin[y]).append(bin_idx_regions[i,1]))
                else:
                    label_line_bin[y] = tuple(list(label_line_bin[y]).append(bin_regions[i,2]))
        ax.plot(*zip(*plot_points[s]), c=color, marker=".", markersize=10-(n_bins / 10), zorder=zorder)
    
    # Plot axes
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    if subplot:
        if args.evenly_spaced_intervals:
            ax2.set_xlabel("Interval index")
        else:
            ax2.set_xlabel("Genomic position")
    else:
        if args.evenly_spaced_intervals:
            ax.set_xlabel("Interval index")
        else:
            ax.set_xlabel("Genomic position")
    ax.set_ylabel("Copy ratio")
    label_x = x_max + x_len * 0.025
    
    # Display sample names
    if args.show_samples:
        sample_label_y_buffer = (y_max - y_min) * 0.05
        if subplot:
            sample_label_y_buffer += sample_label_y_buffer * (1 / subplot_ratio)
        label_clusters = cluster_labels(sample_labels, sample_label_y_buffer)
        for c_min, c_max in label_clusters:
            old_y_list = sorted(list(label_clusters[(c_min, c_max)]))
            for j in range(len(old_y_list)):
                ax.text(label_x, c_min + sample_label_y_buffer * (j + 0.5), sample_labels[old_y_list[j]])
                line_y = [old_y_list[j], c_min + sample_label_y_buffer * (j + 0.7)]
                for x in label_line_bin[old_y_list[j]]:
                    line_x = [x + x_buffer * 0.375, label_x - x_buffer * 0.125]
                    ax.plot(line_x, line_y, c="0", clip_on=False, linewidth=1)
    
    # Display functional annotations 
    if subplot:
        ax2.set_ylim(0, 1)
        label_heights = np.arange(1, len(subplot_labels) + 1) / (len(subplot_labels) + 1) - (0.2 / (len(subplot_labels) + 1))
        if len(subplot_labels) > 1:
            patch_height = 0.4 / (len(subplot_labels) + 1)
        else:
            patch_height = 0.25
        if args.max_annotation_intensities is None:
            cbar_array = [np.arange(25, 176) / 200]
        else:
            cbar_array_interval = min(args.max_annotation_intensities, 151)
            cbar_array = [[(25 + i / (cbar_array_interval - 1) * 150) / 200 for i in range(cbar_array_interval)]]
        for i in range(len(subplot_labels)):
            ax2.text(label_x, label_heights[i], subplot_labels[i])
            cbar_left = label_x + 0.325 * x_len
            cbar_right = label_x + 0.425 * x_len
            cbar_bottom = label_heights[i]
            cbar_top = cbar_bottom + patch_height
            plt.imshow(cbar_array, cmap="Greys", vmin=0, vmax=1, aspect="auto", extent=(cbar_left, cbar_right, cbar_bottom, cbar_top), clip_on=False, zorder=0)
            if subplot_labels[i] == "Missing fraction":
                mat = np.mean(np.isnan(cr_matrix), axis=1)
                max_missing = np.nanmax(mat)
                min_missing = np.nanmin(mat)
                range_missing = max_missing - min_missing
                ax2.text(label_x + 0.2625 * x_len, cbar_bottom, str("{:.2f}".format(min_missing)), size=8, zorder=1)
                ax2.text(label_x + 0.4375 * x_len, cbar_bottom, str("{:.2f}".format(max_missing)), size=8, zorder=1)
                for r in range(n_regions):
                    if range_missing > 0:
                        color = str(0.75 * (1 - ((mat[r] - min_missing) / range_missing)) + 0.125)
                        if args.max_annotation_intensities is not None:
                            for j in range(args.max_annotation_intensities):
                                cint_start = min_missing + range_missing * j / (args.max_annotation_intensities)
                                cint_end = min_missing + range_missing * (j + 1) / (args.max_annotation_intensities)
                                if float(color) >= cint_start and float(color) <= cint_end:
                                    color = str((0.125 + j / (args.max_annotation_intensities - 1) * 0.75))
                                    break
                    else:
                        color = "0.625"
                    if args.evenly_spaced_intervals:
                        for j in range(n_bins):
                            if r >= bin_idx_regions[j,0] and r < bin_idx_regions[j,1]:
                                jbin_len = bin_regions[j,2] - bin_regions[j,1]
                                jbin_regs = bin_idx_regions[j,1] - bin_idx_regions[j,0]
                                left_coord = (max(bin_regions[j,1], cr_regions[r,1]) - bin_regions[j,1]) / jbin_len * jbin_regs + bin_idx_regions[j,0]
                                right_coord = (min(bin_regions[j,2], cr_regions[r,2]) - bin_regions[j,1]) / jbin_len * jbin_regs + bin_idx_regions[j,0]
                                ax2.add_patch(Rectangle((left_coord, label_heights[i]), right_coord - left_coord, patch_height, color=color, zorder=0))
                                break
                    else:
                        ax2.add_patch(Rectangle((cr_regions[r,1], label_heights[i]), cr_regions[r,2] - cr_regions[r,1], patch_height, color=color, zorder=0))
            else:
                if subplot_labels[i] == "GC content":
                    d = gc_dict
                elif subplot_labels[i] == "Segdup content":
                    d = sd_dict
                elif subplot_labels[i] == "Mappability":
                    d = m_dict
                max_d = max(d.values())
                min_d = min(d.values())
                range_d = max_d - min_d
                ax2.text(label_x + 0.2625 * x_len, cbar_bottom, str("{:.2f}".format(min_d)), size=8, zorder=1)
                ax2.text(label_x + 0.4375 * x_len, cbar_bottom, str("{:.2f}".format(max_d)), size=8, zorder=1)
                for s, e in d:
                    if range_d > 0:
                        color = str(0.75 * (1 - ((d[(s,e)] - min_d) / range_d)) + 0.125)
                        if args.max_annotation_intensities is not None:
                            for j in range(args.max_annotation_intensities):
                                cint_start = min_d + range_d * j / (args.max_annotation_intensities)
                                cint_end = min_d + range_d * (j + 1) / (args.max_annotation_intensities)
                                if float(color) >= cint_start and float(color) <= cint_end:
                                    color = str((0.125 + j / (args.max_annotation_intensities - 1) * 0.75))
                                    break
                    else:
                        color = "0.625"
                    if args.evenly_spaced_intervals:
                        for j in range(n_bins):
                            if e >= bin_regions[j,1] and s <= bin_regions[j,2]:
                                jbin_len = bin_regions[j,2] - bin_regions[j,1]
                                jbin_regs = bin_idx_regions[j,1] - bin_idx_regions[j,0]
                                left_coord = (max(bin_regions[j,1], s) - bin_regions[j,1]) / jbin_len * jbin_regs + bin_idx_regions[j,0]
                                right_coord = (min(bin_regions[j,2], e) - bin_regions[j,1]) / jbin_len * jbin_regs + bin_idx_regions[j,0]
                                ax2.add_patch(Rectangle((left_coord, label_heights[i]), right_coord - left_coord, patch_height, color=color, zorder=0))
                    else:
                        ax2.add_patch(Rectangle((s, label_heights[i]), e - s, patch_height, color=color, zorder=0))
    
    # Output plot
    plt.show()
    fig.savefig(args.output, bbox_inches="tight")
    
if __name__ == "__main__":
    
    import sys
    import argparse
    
    # Inputs
    parser = argparse.ArgumentParser(description="Plot copy ratios from whole exome sequencing data.")
    parser.add_argument("data", help="path to a regions x samples tabix matrix containing copy ratio data")
    parser.add_argument("chrom", help="the chromosome to plot")
    parser.add_argument("--start", "-s", default=0, type=int, help="the genomic coordinate at which to start plotting")
    parser.add_argument("--end", "-e", default=float("inf"), type=int, help="the genomic coordinate at which to finish plotting")
    parser.add_argument("--auto_highlight", "-ah", action="store_true", help="automatically highlight samples with directionally consistent aberrant copy ratios in at least 2 adjacent bins")
    parser.add_argument("--exp_cr", "-exp", default=2, type=float, help="the expected copy ratio")
    parser.add_argument("--aberrant_threshold", "-t", default=1, type=float, help="the deviation threshold at which copy ratios are deemed aberrant")
    parser.add_argument("--aberrant_dup_color", "-dupc", default="b", help="the color to highlight aberrantly high copy ratios")
    parser.add_argument("--aberrant_del_color", "-delc", default="r", help="the color to highlight aberrantly low copy ratios")
    parser.add_argument("--show_samples", "-S", action="store_true", help="add the names of highlighted samples on the plot")
    parser.add_argument("--sample_highlight", "-sh", nargs="*", default=[], help="specific sample names to highlight")
    parser.add_argument("--sample_color", "-sc", default="indigo", help="the color to highlight specified sample names")
    parser.add_argument("--show_intervals", "-i", action="store_true", help="shade intervals corresponding to exons on the main plot")
    parser.add_argument("--show_missing", "-M", action="store_true", help="display the fraction of missing samples in each interval")
    parser.add_argument("--show_gc", "-gc", action="store_true", help="display the GC content of intervals")
    parser.add_argument("--show_segdup", "-sd", action="store_true", help="display the segdup content of intervals")
    parser.add_argument("--show_map", "-m", action="store_true", help="display the mappability of intervals")
    parser.add_argument("--functional_annotations", "-f", help="a .tsv file containing the GC content, segdup content, and mappability of each interval (required if --show_gc, --show_segdup, or --show_map is enabled)")
    parser.add_argument("--max_annotation_intensities", "-ai", type=int, help="the maximum number of different intensities shown for annotations")
    parser.add_argument("--evenly_spaced_intervals", "-I", action="store_true", help="display intervals as adjacent equal size regions rather than by genomic coordinates")
    parser.add_argument("--max_bins", "-b", default=float("inf"), type=int, help="the maximum number of interval bins to show in the plot; if intervals must be pooled into bins, the mean copy ratio for each sample will be displayed for each bin")
    parser.add_argument("--evenly_pooled_bins", "-P", action="store_true", help="define bins to contain equal numbers of intervals rather than equal genomic lengths")
    parser.add_argument("--bin_median", "-bm", action="store_true", help="plot the median copy ratio for each bin rather than the mean copy ratio")
    parser.add_argument("--output", "-o", default="copy_ratio_plot.png", help="output file where the resulting plot will be saved")
    parser.add_argument("--last_path", "-lp", help="moves a path to the end of sys.path; for use in virtual environments on ERISOne")
    args = parser.parse_args()
    
    # Path modifications for virtual environments on ERISOne
    if args.last_path is not None and args.last_path in sys.path:
        sys.path.remove(args.last_path)
        sys.path.append(args.last_path)
    
    import pysam
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib.patches import Rectangle
    
    create_plot(args)
    
    