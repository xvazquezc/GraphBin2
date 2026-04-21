#!/usr/bin/env python3

import concurrent.futures
import csv
import heapq
import itertools as it
import logging
import os
import re
import subprocess
import sys
import time

import multiprocessing
from collections import deque

from cogent3.parse.fasta import MinimalFastaParser
from igraph import *
from tqdm import tqdm

from .bidirectionalmap.bidirectionalmap import BidirectionalMap

__author__ = "Vijini Mallawaarachchi, Anuradha Wickramarachchi, and Yu Lin"
__copyright__ = "Copyright 2020, GraphBin2 Project"
__license__ = "BSD"
__version__ = "1.3.3"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Stable Release"

# Shared state for ProcessPoolExecutor workers (fork-inherited on Linux).
# Set in run() before the pool is created.
_is_multi_kwargs: dict = {}


def _is_multi_worker(n: int):
    """Module-level worker callable for ProcessPoolExecutor."""
    return is_multi(contig=n, **_is_multi_kwargs)


def run(args):
    # Get arguments
    # ---------------------------------------------------

    contigs_file = args.contigs
    assembly_graph_file = args.graph
    contig_paths = args.paths
    contig_bins_file = args.binned
    output_path = args.output
    prefix = args.prefix
    depth = args.depth
    threshold = args.threshold
    delimiter = args.delimiter
    nthreads = args.nthreads

    n_bins = 0

    # Setup logger
    # -----------------------
    logger = logging.getLogger(f"GraphBin2 {__version__}")
    logger.setLevel(logging.DEBUG)
    if not logger.handlers:
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        consoleHeader = logging.StreamHandler()
        consoleHeader.setFormatter(formatter)
        consoleHeader.setLevel(logging.INFO)
        logger.addHandler(consoleHeader)

    # Setup output path for log file
    # ---------------------------------------------------

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    fileHandler = logging.FileHandler(f"{output_path}/{prefix}graphbin2.log")
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    logger.info(
        "Welcome to GraphBin2: Refined and Overlapped Binning of Metagenomic Contigs using Assembly Graphs."
    )
    logger.info(
        "This version of GraphBin2 makes use of the assembly graph produced by MEGAHIT which is based on the de Bruijn graph approach."
    )

    logger.info(f"Input arguments:")
    logger.info(f"Contigs file: {contigs_file}")
    logger.info(f"Assembly graph file: {assembly_graph_file}")
    logger.info(f"Contig paths file: {contig_paths}")
    logger.info(f"Existing binning output file: {contig_bins_file}")
    logger.info(f"Final binning output file: {output_path}")
    logger.info(f"Depth: {depth}")
    logger.info(f"Threshold: {threshold}")
    logger.info(f"Number of threads: {nthreads}")

    logger.info(f"GraphBin2 started")

    start_time = time.time()

    # Get length and coverage of contigs
    # --------------------------------------------------------

    original_contigs = {}
    contig_descriptions = {}

    for label, seq in MinimalFastaParser(contigs_file):
        name = label.split()[0]
        original_contigs[name] = seq
        contig_descriptions[name] = label

    # Build assembly graph
    # -------------------------------------

    node_count = 0

    graph_contigs = {}

    links = []

    my_map = BidirectionalMap()

    try:
        # Get links from .gfa file
        with open(assembly_graph_file) as file:
            for line in file.readlines():
                line = line.strip()

                # Identify lines with link information
                if line.startswith("L"):
                    link = []

                    strings = line.split("\t")

                    start_1 = "NODE_"
                    end_1 = "_length"

                    link1 = int(
                        re.search("%s(.*)%s" % (start_1, end_1), strings[1]).group(1)
                    )

                    start_2 = "NODE_"
                    end_2 = "_length"

                    link2 = int(
                        re.search("%s(.*)%s" % (start_2, end_2), strings[3]).group(1)
                    )

                    link.append(link1)
                    link.append(link2)
                    links.append(link)

                elif line.startswith("S"):
                    strings = line.split()

                    start = "NODE_"
                    end = "_length"

                    contig_num = int(
                        re.search("%s(.*)%s" % (start, end), strings[1]).group(1)
                    )

                    my_map[node_count] = int(contig_num)

                    graph_contigs[contig_num] = strings[2]

                    node_count += 1

        logger.info(f"Total number of contigs available: {node_count}")

        contigs_map = my_map
        contigs_map_rev = my_map.inverse

        # Create graph
        assembly_graph = Graph()

        # Add vertices
        assembly_graph.add_vertices(node_count)

        # Create list of edges
        edge_list = []

        for i in range(node_count):
            assembly_graph.vs[i]["id"] = i
            assembly_graph.vs[i]["label"] = str(contigs_map[i])

        # Iterate links
        for link in links:
            # Remove self loops
            if link[0] != link[1]:
                # Add edge to list of edges
                edge_list.append((contigs_map_rev[link[0]], contigs_map_rev[link[1]]))

        # Add edges to the graph
        assembly_graph.add_edges(edge_list)
        assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            f"Please make sure that the correct path to the assembly graph file is provided."
        )
        logger.info(f"Exiting GraphBin2... Bye...!")
        sys.exit(1)

    logger.info(f"Total number of edges in the assembly graph: {len(edge_list)}")

    # Map original contig IDs to contig IDS of assembly graph
    # --------------------------------------------------------

    graph_to_contig_map = BidirectionalMap()

    # Build reverse lookup: sequence -> original contig name
    original_seq_to_name = {}
    for name, seq in original_contigs.items():
        original_seq_to_name[seq] = name

    for graph_name, graph_seq in graph_contigs.items():
        if graph_seq in original_seq_to_name:
            graph_to_contig_map[graph_name] = original_seq_to_name[graph_seq]

    graph_to_contig_map_rev = graph_to_contig_map.inverse

    # Get length and coverage of contigs
    # --------------------------------------------------------

    contig_lengths = {}
    coverages = {}

    my_map = BidirectionalMap()

    for label, seq in MinimalFastaParser(contigs_file):
        strings = label.split()

        contig_num = contigs_map_rev[graph_to_contig_map_rev[strings[0]]]
        length = int(strings[3][4:])
        coverage = int(float(strings[2][6:]))

        contig_lengths[contig_num] = length
        coverages[contig_num] = coverage

    # Get the number of bins from the initial binning result
    # --------------------------------------------------------

    try:
        all_bins_list = []

        with open(contig_bins_file) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=delimiter)
            for row in readCSV:
                all_bins_list.append(row[1])

        bins_list = list(set(all_bins_list))
        bins_list.sort()

        n_bins = len(bins_list)
        logger.info(f"Number of bins available in binning result: {n_bins}")

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            f"Please make sure that the correct path to the binning result file is provided and it is having the correct format"
        )
        logger.info(f"Exiting GraphBin2... Bye...!")
        sys.exit(1)

    # Get initial binning result
    # ----------------------------

    bins = [set() for x in range(n_bins)]

    try:
        with open(contig_bins_file) as contig_bins:
            readCSV = csv.reader(contig_bins, delimiter=delimiter)
            for row in readCSV:
                contig_num = contigs_map_rev[int(graph_to_contig_map_rev[row[0].split()[0]])]

                bin_num = bins_list.index(row[1])
                bins[bin_num].add(contig_num)

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            f"Please make sure that you have provided the correct assembler type and the correct path to the binning result file in the correct format."
        )
        logger.info(f"Exiting GraphBin2... Bye...!")
        sys.exit(1)

    # Get binned and unbinned contigs
    # -----------------------------------------------------

    binned_contigs = set()
    for n in range(n_bins):
        binned_contigs.update(bins[n])

    unbinned_contigs = set(range(node_count)) - binned_contigs

    # Reverse map: contig_id → bin index for O(1) lookup throughout
    contig_bin_map = {c: n for n in range(n_bins) for c in bins[n]}

    logger.info(f"Number of binned contigs: {len(binned_contigs)}")
    logger.info(f"Total number of unbinned contigs: {len(unbinned_contigs)}")

    # Get isolated vertices
    # -----------------------------------------------------

    isolated = set()

    for i in range(node_count):
        neighbours = assembly_graph.neighbors(i, mode=ALL)

        if len(neighbours) == 0:
            isolated.add(i)

    logger.info(f"Number of isolated contigs: {len(isolated)}")

    # The BFS function to search labelled nodes
    # -----------------------------------------------------

    def runBFS(node, threhold=depth):
        queue = deque([node])
        visited = set()
        depth = {}

        depth[node] = 0

        labelled_nodes = set()

        while queue:
            active_node = queue.popleft()
            visited.add(active_node)

            if active_node in binned_contigs and len(visited) > 1:
                contig_bin = contig_bin_map.get(active_node, -1)
                if contig_bin != -1:
                    labelled_nodes.add(
                        (
                            node,
                            active_node,
                            contig_bin,
                            depth[active_node],
                            abs(
                                coverages.get(node, 0) - coverages.get(active_node, 0)
                            ),
                        )
                    )

            else:
                for neighbour in assembly_graph.neighbors(active_node, mode=ALL):
                    if neighbour not in visited:
                        depth[neighbour] = depth[active_node] + 1
                        if depth[neighbour] > threhold:
                            continue
                        queue.append(neighbour)

        return labelled_nodes

    # Remove labels of unsupported vertices
    # -----------------------------------------------------

    logger.info(f"Removing labels of unsupported vertices")

    iter_num = 1

    while True:
        logger.debug(f"Iteration: {iter_num}")

        remove_labels = {}

        # Initialise progress bar
        pbar = tqdm(total=len(binned_contigs))

        for my_node in binned_contigs:
            if my_node not in isolated:
                my_contig_bin = contig_bin_map.get(my_node, -1)

                BFS_labelled_nodes = list(runBFS(my_node))

                if len(BFS_labelled_nodes) > 0:
                    # Get the count of nodes in the closest_neighbours that belongs to each bin
                    BFS_labelled_bin_counts = [0 for x in range(n_bins)]

                    for i in range(len(BFS_labelled_nodes)):
                        BFS_labelled_bin_counts[BFS_labelled_nodes[i][2]] += 1

                    zero_bin_count = 0

                    # Count the number of bins which have no BFS_labelled_contigs
                    for j in BFS_labelled_bin_counts:
                        if j == 0:
                            zero_bin_count += 1

                    # Get the bin number which contains the maximum number of BFS_labelled_contigs
                    max_index = BFS_labelled_bin_counts.index(
                        max(BFS_labelled_bin_counts)
                    )

                    # If there are no BFS nodes of same label as contig, remove label
                    if (
                        my_contig_bin != -1
                        and BFS_labelled_bin_counts[my_contig_bin] == 0
                    ):
                        remove_labels[my_node] = my_contig_bin

                    # Check if all the BFS_labelled_contigs are in one bin
                    elif zero_bin_count == (len(BFS_labelled_bin_counts) - 1):
                        # If contig is not in the bin with maximum number of BFS_labelled_contigs
                        if (
                            max_index != my_contig_bin
                            and BFS_labelled_bin_counts[max_index] > 1
                            and contig_lengths.get(my_node, 10000) < 10000
                        ):
                            remove_labels[my_node] = my_contig_bin

            # Update progress bar
            pbar.update(1)

        # Close progress bar
        pbar.close()

        if len(remove_labels) == 0:
            break
        else:
            for contig in remove_labels:
                bins[remove_labels[contig]].discard(contig)
                binned_contigs.discard(contig)
                unbinned_contigs.add(contig)
                contig_bin_map.pop(contig, None)

        iter_num += 1

    # Refine labels of inconsistent vertices
    # -----------------------------------------------------

    logger.info(f"Refining labels of inconsistent vertices")

    iter_num = 1

    once_moved = set()

    while True:
        logger.debug(f"Iteration: {iter_num}")

        contigs_to_correct = {}

        # Initialise progress bar
        pbar = tqdm(total=len(binned_contigs))

        for my_node in binned_contigs:
            if my_node not in isolated and my_node not in once_moved:
                my_contig_bin = contig_bin_map.get(my_node, -1)

                BFS_labelled_nodes = list(runBFS(my_node))

                # Get the count of nodes in the closest_neighbours that belongs to each bin
                BFS_labelled_bin_counts = [0 for x in range(n_bins)]

                for i in range(len(BFS_labelled_nodes)):
                    BFS_labelled_bin_counts[BFS_labelled_nodes[i][2]] += 1

                zero_bin_count = 0

                # Count the number of bins which have no BFS_labelled_contigs
                for j in BFS_labelled_bin_counts:
                    if j == 0:
                        zero_bin_count += 1

                # Get the bin number which contains the maximum number of BFS_labelled_contigs
                max_index = BFS_labelled_bin_counts.index(max(BFS_labelled_bin_counts))

                weighted_bin_count = [0 for x in range(n_bins)]

                for i in range(len(BFS_labelled_nodes)):
                    path_length = BFS_labelled_nodes[i][3]
                    weighted_bin_count[BFS_labelled_nodes[i][2]] += 1 / (2**path_length)

                should_move = False

                max_weight = -1
                max_weight_bin = -1

                for i in range(len(weighted_bin_count)):
                    if (
                        len(BFS_labelled_nodes) > 0
                        and my_contig_bin != -1
                        and i != my_contig_bin
                        and weighted_bin_count[i] > 0
                        and weighted_bin_count[i]
                        > weighted_bin_count[my_contig_bin] * threshold
                    ):
                        should_move = True
                        if max_weight < weighted_bin_count[i]:
                            max_weight = weighted_bin_count[i]
                            max_weight_bin = i

                if should_move and max_weight_bin != -1:
                    contigs_to_correct[my_node] = (my_contig_bin, max_weight_bin)
                    once_moved.add(my_node)

            # Update progress bar
            pbar.update(1)

        # Close progress bar
        pbar.close()

        if len(contigs_to_correct) == 0:
            break
        else:
            for contig in contigs_to_correct:
                old_bin = contigs_to_correct[contig][0]
                new_bin = contigs_to_correct[contig][1]
                bins[old_bin].discard(contig)
                bins[new_bin].add(contig)
                contig_bin_map[contig] = new_bin

        iter_num += 1

    # Get non isolated contigs

    logger.info(f"Obtaining non isolated contigs")

    # Initialise progress bar
    pbar = tqdm(total=node_count)

    non_isolated_set = set()
    for component in assembly_graph.clusters():
        if any(m in binned_contigs for m in component):
            non_isolated_set.update(component)
        pbar.update(len(component))
    pbar.close()

    non_isolated = list(non_isolated_set)

    logger.info(f"Number of non-isolated contigs: {len(non_isolated)}")

    non_isolated_unbinned = non_isolated_set.intersection(unbinned_contigs)

    logger.info(
        f"Number of non-isolated unbinned contigs: {len(non_isolated_unbinned)}"
    )

    # Propagate labels to unlabelled vertices
    # -----------------------------------------------------

    logger.info(f"Propagating labels to unlabelled vertices")

    # Initialise progress bar
    pbar = tqdm(total=len(non_isolated_unbinned))

    contigs_to_bin = set()

    for contig in binned_contigs:
        if contig in non_isolated_set:
            closest_neighbours = filter(
                lambda x: x not in binned_contigs,
                assembly_graph.neighbors(contig, mode=ALL),
            )
            contigs_to_bin.update(closest_neighbours)

    # Improvement 4: heap entries are plain (dist, cov_diff, to_bin, src, bin_) tuples;
    # Python's native tuple comparison replaces the DataWrap wrapper class.
    sorted_node_list = []
    for x in contigs_to_bin:
        for nd, src, bin__, bfs_dist, cov_diff in runBFS(x, threhold=depth):
            heapq.heappush(sorted_node_list, (bfs_dist, cov_diff, nd, src, bin__))

    # Improvement 1: lazy-deletion drain — the filter+heapify block is removed.
    # Stale entries (already-assigned nodes) are discarded at pop-time by the
    # `if to_bin in non_isolated_unbinned` set-membership check (O(1)).
    #
    # Improvement 2: when a node is assigned, its direct unbinned neighbours get a
    # single distance-1 candidate pushed directly from the just-assigned node, instead
    # of a full runBFS call. Any path to those neighbours from previously-labelled nodes
    # is already in the heap, so this push supplies the only new information.
    while sorted_node_list:
        bfs_dist, cov_diff, to_bin, src, bin_ = heapq.heappop(sorted_node_list)

        if to_bin in non_isolated_unbinned:
            bins[bin_].add(to_bin)
            binned_contigs.add(to_bin)
            non_isolated_unbinned.discard(to_bin)
            unbinned_contigs.discard(to_bin)
            contig_bin_map[to_bin] = bin_

            # Update progress bar
            pbar.update(1)

            # Push distance-1 candidates for newly reachable unbinned neighbours.
            for n in assembly_graph.neighbors(to_bin, mode=ALL):
                if n not in binned_contigs:
                    cov_d = abs(coverages.get(n, 0) - coverages.get(to_bin, 0))
                    heapq.heappush(sorted_node_list, (1, cov_d, n, to_bin, bin_))

    # Close progress bar
    pbar.close()

    # Determine contigs belonging to multiple bins
    # -----------------------------------------------------

    logger.info(f"Determining multi-binned contigs")

    bin_cov_sum = [0 for x in range(n_bins)]
    bin_contig_len_total = [0 for x in range(n_bins)]

    for i in range(n_bins):
        for contig in bins[i]:
            if contig in non_isolated_set:
                node_cov = coverages.get(contig, 0)
                node_len = contig_lengths.get(contig, 0)
                bin_cov_sum[i] += node_cov * node_len
                bin_contig_len_total[i] += node_len

    # Set up shared state for worker processes (fork-inherited on Linux)
    global _is_multi_kwargs
    _is_multi_kwargs = dict(
        non_isolated=non_isolated_set,
        binned_contigs=binned_contigs,
        n_bins=n_bins,
        bins=bins,
        bin_cov_sum=bin_cov_sum,
        bin_contig_len_total=bin_contig_len_total,
        coverages=coverages,
        contigs_map=contigs_map,
        contig_lengths=contig_lengths,
        assembly_graph=assembly_graph,
        contig_bin_map=contig_bin_map,
    )
    chunksize = max(1, node_count // (nthreads * 4))
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=nthreads,
        mp_context=multiprocessing.get_context("fork"),
    ) as executor:
        mapped = list(
            tqdm(
                executor.map(_is_multi_worker, range(node_count), chunksize=chunksize),
                total=node_count,
            )
        )

    # Get multi-bin results
    multi_bins = list(filter(lambda x: x is not None, mapped))

    if len(multi_bins) == 0:
        logger.info(f"No multi-labelled contigs were found ==>")
    else:
        logger.info(f"Found {len(multi_bins)} multi-labelled contigs")

    # Add contigs to multiplt bins
    for contig, min_diff_combination in multi_bins:
        contig_label = graph_to_contig_map.get(contigs_map[contig], str(contigs_map[contig]))
        logger.info(
            contig_label
            + " belongs to bins "
            + ", ".join(bins_list[s] for s in min_diff_combination)
        )
        for mybin in min_diff_combination:
            bins[mybin].add(contig)

    # Determine elapsed time
    elapsed_time = time.time() - start_time

    # Show elapsed time for the process
    logger.info(f"Elapsed time: {elapsed_time} seconds")

    # Write result to output file
    # -----------------------------------

    logger.info(f"Writing the final binning results to file")

    output_bins = []

    final_bins = {}

    for i in range(n_bins):
        for contig in bins[i]:
            final_bins[contig] = bins_list[i]

    output_bins_path = f"{output_path}{prefix}bins/"

    if not os.path.isdir(output_bins_path):
        subprocess.run(f"mkdir -p {output_bins_path}", shell=True)

    bin_files = {}

    for bin_num in range(n_bins):
        bin_files[bins_list[bin_num]] = open(
            f"{output_bins_path}{prefix}{bins_list[bin_num]}.fasta", "w+"
        )

    for label, seq in MinimalFastaParser(contigs_file):
        strings = label.split()
        orig_name = strings[0]
        if orig_name not in graph_to_contig_map_rev:
            continue
        contig_num = contigs_map_rev[graph_to_contig_map_rev[orig_name]]

        if contig_num in final_bins:
            bin_files[final_bins[contig_num]].write(f">{label}\n{seq}\n")

    # Close output files
    for bin_num in range(n_bins):
        bin_files[bins_list[bin_num]].close()

    for k in range(n_bins):
        for contig in bins[k]:
            contig_label = graph_to_contig_map.get(contigs_map[contig])
            if contig_label is None:
                continue
            output_bins.append([contig_label, bins_list[k]])

    output_file = f"{output_path}{prefix}graphbin2_output.csv"

    with open(output_file, mode="w") as output_file:
        output_writer = csv.writer(
            output_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL
        )

        for row in output_bins:
            output_writer.writerow(row)

    logger.info(f"Final binning results can be found at {output_file.name}")

    # Exit program
    # -----------------------------------

    logger.info(f"Thank you for using GraphBin2!")


def is_multi(
    contig,
    non_isolated,
    binned_contigs,
    n_bins,
    bins,
    bin_cov_sum,
    bin_contig_len_total,
    coverages,
    contigs_map,
    contig_lengths,
    assembly_graph,
    contig_bin_map=None,
):
    if contig in non_isolated and contig in binned_contigs:
        # O(1) bin lookup when contig_bin_map is available
        if contig_bin_map is not None:
            contig_bin = contig_bin_map.get(contig, -1)
        else:
            contig_bin = next(
                (n for n in range(n_bins) if contig in bins[n]), -1
            )

        # Get average coverage of each connected component representing a bin excluding the contig
        bin_coverages = list(bin_cov_sum)
        bin_contig_lengths = list(bin_contig_len_total)

        bin_coverages[contig_bin] = bin_coverages[contig_bin] - (
            coverages.get(contig, 0) * contig_lengths.get(contig, 0)
        )
        bin_contig_lengths[contig_bin] = (
            bin_contig_lengths[contig_bin] - contig_lengths.get(contig, 0)
        )

        for i in range(n_bins):
            if bin_contig_lengths[i] != 0:
                bin_coverages[i] = bin_coverages[i] / bin_contig_lengths[i]

        # Get coverages of neighbours
        neighbour_bins = [[] for x in range(n_bins)]

        neighbour_bin_coverages = [[] for x in range(n_bins)]

        neighbours = assembly_graph.neighbors(contig, mode=ALL)

        for neighbour in neighbours:
            # O(1) neighbour bin lookup
            if contig_bin_map is not None:
                nb = contig_bin_map.get(neighbour, -1)
                if nb != -1:
                    neighbour_bins[nb].append(neighbour)
                    neighbour_bin_coverages[nb].append(coverages.get(neighbour, 0))
            else:
                for n in range(n_bins):
                    if neighbour in bins[n]:
                        neighbour_bins[n].append(neighbour)
                        neighbour_bin_coverages[n].append(coverages.get(neighbour, 0))
                        break

        zero_bin_count = 0

        non_zero_bins = []

        # Count the number of bins which have no labelled neighbouring contigs
        for i in range(len(neighbour_bins)):
            if len(neighbour_bins[i]) == 0:
                zero_bin_count += 1
            else:
                non_zero_bins.append(i)

        if zero_bin_count < n_bins - 1:
            bin_combinations = []

            for i in range(len(non_zero_bins)):
                bin_combinations += list(it.combinations(non_zero_bins, i + 1))

            min_diff = sys.maxsize
            min_diff_combination = -1

            for combination in bin_combinations:
                comb_cov_total = 0

                for i in range(len(combination)):
                    comb_cov_total += bin_coverages[combination[i]]

                cov_diff = abs(comb_cov_total - coverages.get(contig, 0))

                if cov_diff < min_diff:
                    min_diff = cov_diff
                    min_diff_combination = combination

            if (
                min_diff_combination != -1
                and len(min_diff_combination) > 1
                and contig_lengths.get(contig, 0) > 1000
            ):
                # return True
                return contig, min_diff_combination

    return None


def main(args):
    run(args)


if __name__ == "__main__":
    main()
