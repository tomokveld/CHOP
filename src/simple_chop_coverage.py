from collections import defaultdict
from intervaltree import IntervalTree
import argparse
import itertools
# import matplotlib.pyplot as plt
import numpy as np
import pysam
import re
import sys


# Very simple example how to project a pileup of coverage back onto the graph (alternate alignments are not included here!)
# Input indexed BAM files -> samtools view -bS input.sam | samtools sort -T TMP_FILE -o input.bam ; samtools index input.bam

class DefaultHelpParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        sys.exit(2)


def gen_interval_tree(interval):
    # Returns an interval tree queryable by points and ranges
    # [0, n] -> (node, offset)
    #     Node        Node        Node
    # |---------|--------------|----------|

    offset = 0

    tree = IntervalTree()
    for pair in interval.split('|', 1)[-1].split('_'):
        (node, start, end) = map(int, re.split('[+,]', pair))

        tree[start + offset:end + offset] = (node, offset)
        offset += end

    return (tree, tree.begin())


def prev_cur_and_next(iterable):
    previous, items, following = itertools.tee(iterable, 3)
    previous = itertools.chain([None], previous)
    following = itertools.chain(itertools.islice(following, 1, None), [None])
    return itertools.izip(previous, items, following)


def is_array_gapped(xs):
    for p, c, n in prev_cur_and_next(xs):
        if not n:
            if c - p != 1:
                return True
        else:
            if n - c != 1:
                return True
    return False


# Nicer to do this without pre-initialization
def init_coverage_dict(in_graph):
    coverage_dict = defaultdict(list)

    with open(in_graph, 'rb') as f_in:
        for line in f_in:
            if line.startswith('S'):
                (_, node, seq) = line.split()[:3]
                coverage_dict[int(node)] = np.zeros(len(seq))

    return coverage_dict


def init_remap_dict(in_map):
    if in_map:
        interval_dict = defaultdict(str)

        with open(in_map, 'rb') as f_in:
            for line in f_in:
                key, interval = line.split()
                interval_dict[key] = interval.split('|', 1)[-1]
        return interval_dict
    else:
        return 0


def gen_alignment_coverage(in_aln, coverage_dict, interval_dict):
    with pysam.AlignmentFile(in_aln, 'rb') as aln_in:
        try:
            aln_in.check_index()
        except ValueError:
            print "No BAM index found!"
            exit(2)

        for interval, length in zip(aln_in.references, aln_in.lengths):
            if interval_dict:
                (tree, offset) = gen_interval_tree(interval_dict[interval])
            else:
                (tree, offset) = gen_interval_tree(interval)

            # xs = []

            prev_node = None
            prev_nsegments = None
            prev_pileuppos = None

            # Extract read pileup over full length of interval
            for pileup in aln_in.pileup(interval, 0, length):
                tree_offset = offset + pileup.pos
                qs = tree[tree_offset].pop()

                node = qs.data[0]
                node_offset = tree_offset - qs.data[1]

                if node == 1:
                    continue

                # Assign node coverage
                coverage_dict[node][node_offset] = pileup.nsegments

                # Assign edge coverage
                if prev_node:
                    if prev_node != node and pileup.pos - prev_pileuppos == 1:
                        if not coverage_dict[prev_node, node]:
                            coverage_dict[prev_node, node] = min(
                                prev_nsegments, pileup.nsegments)
                        else:
                            coverage_dict[prev_node,
                                          node] += min(prev_nsegments, pileup.nsegments)

                prev_node = node
                prev_nsegments = pileup.nsegments
                prev_pileuppos = pileup.pos

                # xs.append(pileup.pos)

            # print is_array_gapped(xs)
            # plt.scatter(xs, xs)
            # plt.show()
            # plt.clf()

    return coverage_dict


def write_coverage(coverage_dict, summation):
    for key, value in coverage_dict.iteritems():
        if isinstance(value, int):
            print "{}\t{}".format(key, value)
        else:
            if summation:
                print "{}\t{}".format(key, sum(value) / len(value))
            else:
                print "{}\t{}".format(key, value.tolist())


if __name__ == '__main__':
    parser = DefaultHelpParser(
        description="Parse BAM file and output coverage.")
    parser.add_argument('-i', '--bam', dest='in_aln',
                        help='BAM file', type=str, required=True)
    parser.add_argument('-g', '--graph', dest='in_graph',
                        help='GFA graph file', type=str, required=True)
    parser.add_argument('-s', '--sum', dest='sum',
                        help='Sum node coverage (normalized by node length)', required=False, action='store_true')
    parser.add_argument('-r', '--remap', dest='remap',
                        help='Tab seperated file mapping integers to intervals', type=str, default='')

    args = parser.parse_args()

    coverage_dict = gen_alignment_coverage(
        args.in_aln, init_coverage_dict(args.in_graph), init_remap_dict(args.remap))

    write_coverage(coverage_dict, args.sum)
