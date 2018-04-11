from bitarray import bitarray
from collections import defaultdict, deque
from copy import deepcopy
from functools import partial
from itertools import islice, chain, product
from mmap import mmap
from string import maketrans
from time import time
import gzip
import logging
import networkx as nx
import re

logger = logging.getLogger(__name__)

ATTRIBUTE_SET = set(['suffix', 'prefix', 'sequence'])


# TODO: replace lambda with list comprehensions where possible

def profile(f):
    def wrap(*args, **kwargs):
        fname = f.func_name
        argnames = f.func_code.co_varnames[:f.func_code.co_argcount]
        filled_args = ', '.join(
            '%s=%r' % entry
            for entry in zip(argnames, args[:len(argnames)]) + [("args", list(args[len(argnames):]))] + [("kwargs", kwargs)])
        logger.info('Started: {}({})'.format(fname, filled_args))
        starting_time = time()
        output = f(*args, **kwargs)
        logger.info('Ended: Function -> {}, duration {}s'.format(fname,
                                                                 time() - starting_time))
        return output
    return wrap


def skip_header(file_handle, prefix='#'):
    """
    Place the file handle at the position past the header
    """

    last_pos = file_handle.tell()

    while file_handle.readline().startswith(prefix):
        last_pos = file_handle.tell()

    file_handle.seek(last_pos)


def out_deg(graph, nbunch=None):
    if nbunch in graph:
        return len(graph.succ[nbunch])

    if nbunch is None:
        nodes_nbrs = graph.succ.items()
    else:
        nodes_nbrs = ((n, graph.succ[n]) for n in graph.nbunch_iter(nbunch))

    def d_iter():
        for n, nbrs in nodes_nbrs:
            yield (n, len(nbrs))

    return d_iter()


def in_deg(graph, nbunch=None):
    if nbunch in graph:
        return len(graph.pred[nbunch])

    if nbunch is None:
        nodes_nbrs = graph.pred.items()
    else:
        nodes_nbrs = ((n, graph.pred[n]) for n in graph.nbunch_iter(nbunch))

    def d_iter():
        for n, nbrs in nodes_nbrs:
            yield (n, len(nbrs))

    return d_iter()


def convert(text):
    return int(text) if text.isdigit() else text.lower()


def natural_sort(xs):
    return sorted(xs, key=lambda key: [convert(c) for c in re.split('([0-9]+)', key)])


# TODO: Add more file types
def eval_file_type(path):
    magic_dict = {
        "\x1f\x8b\x08": "gz",
        "\x42\x5a\x68": "bz2",
        "\x50\x4b\x03\x04": "zip"
    }

    max_len = max(len(x) for x in magic_dict)

    with open(path) as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    return 0


@profile
def set_to_bit(graph, bits=2):
    haplotypes = sorted(graph.htypes)
    bit_length = len(haplotypes) * bits
    bit_lookup = dict([(j, i * bits) for i, j in enumerate(haplotypes)])
    offset = bits - 1

    for node in graph.nodes():
        idx_ones = set(sum([x() for x in [lambda m=bit_lookup[i]: (
            m, m + offset) for i in graph.node[node]['haplotype']]], ()))
        graph.node[node]['haplotype'] = bitarray(
            ''.join(['1' if i in idx_ones else '0' for i in xrange(bit_length)]))


@profile
def set_to_bit_edges(graph, bits=2):
    haplotypes = sorted(graph.htypes)
    bit_length = len(haplotypes) * bits
    bit_lookup = dict([(j, i * bits) for i, j in enumerate(haplotypes)])
    offset = bits - 1

    for u, v in graph.edges():
        idx_ones = set(sum([x() for x in [lambda m=bit_lookup[i]: (
            m, m + offset) for i in graph[u][v]]], ()))
        graph[u][v] = bitarray(
            ''.join(['1' if i in idx_ones else '0' for i in xrange(bit_length)]))


def count_dfs(graph, u, t):
    if u == t:
        return 1
    else:
        if not graph.node[u]['marked']:
            graph.node[u]['marked'] = sum(
                [count_dfs(graph, c, t) for c in graph.successors(u)])
        return graph.node[u]['marked']


def count_paths(graph):
    # TODO: Do this in a single traversal (start_points and end_points)
    start_points = [node for node, pred in graph.pred.items() if len(
        pred) == 0 and len(graph.successors(node)) > 0]
    end_points = [node for node, succ in graph.succ.items() if len(
        succ) == 0 and len(graph.predecessors(node)) > 0]

    path_len = 0
    for (source, sink) in product(start_points, end_points):
        for node in graph.nodes():
            graph.node[node]['marked'] = 0
        path_len += count_dfs(graph, source, sink)
    return path_len


def eval_dupes(graph):
    identifier_set = set([])
    dupe_dict = defaultdict(lambda: 1)

    for node in graph.nodes():
        ident = generate_identifier(graph.node[node])
        if ident in identifier_set:
            dupe_dict[ident] += 1
        else:
            identifier_set.add(ident)

    return dupe_dict


def cast_attributes(graph):
    for node in graph.nodes():
        for attr in graph.node[node]:
            if not isinstance(graph.node[node][attr], (str, int)):
                graph.node[node][attr] = str(graph.node[node][attr])


def get_n_of_lines(file_path):
    """
    Return total lines in file
    """

    with open(file_path, 'r+') as fp:
        buf = mmap(fp.fileno(), 0)
        lines = 0

        while buf.readline():
            lines += 1

    return lines


def get_node_list(line):
    """
    Reshape list with only specified entries, which are in this case the outgoing and incoming nodes
    """

    # TODO: Use itemgetter instead?
    line = line.split()
    return list(chain([int(line[1])], [int(line[3])]))


def create_connection(graph, line):
    """
    Adds nodes and an edge between them in the graph
    """

    node_list = get_node_list(line)
    graph.add_nodes_from(node_list)
    graph.add_edge(*node_list)


def create_connection_alt(graph, line):
    """
    Adds nodes and an edge between them in the graph
    """

    line = line.split()
    (u, v) = list(chain([int(line[1])], [int(line[3])]))
    graph.add_nodes_from([u, v])
    graph.add_edge(u, v)
    graph[u][v] = set([graph.hdict[i]
                       for i in line[5].split(':')[2].split(';')])


def add_node(graph, line):
    """
    Add a node to the graph with sequence
    """

    line = line.split()
    node = int(line[1])
    seq = line[2]

    if node in graph:
        graph.node[node]['sequence'] = seq
    else:
        graph.add_node(node, sequence=seq)

    # TODO: This assumes that the input we get is always the same...
    if graph.htypes:
        graph.node[node]['haplotype'] = set(
            [graph.hdict[i] for i in line[4].split(':')[2].split(';')])


def add_node_alt(graph, line):
    """
    Add a node to the graph with sequence
    """

    line = line.split()
    node = int(line[1])
    seq = line[2]

    if node in graph:
        graph.node[node]['sequence'] = seq
    else:
        graph.add_node(node, sequence=seq)


def get_header_offset(file_handle, prefix='H'):
    start_pos = file_handle.tell()
    last_pos = start_pos

    while file_handle.readline().startswith(prefix):
        last_pos = file_handle.tell()
    file_handle.seek(start_pos)

    return (start_pos, last_pos)


def get_haplo_header(path):

    with open_gfa(path, 'r') as fin:
        (start, end) = get_header_offset(fin)

        # Not very general
        header = [i for i in (re.split("[\t\:]", fin.read(
            end - start).rstrip().split('\n')[-1])[-1].split(';')) if i]

    return header


# TODO: Generalize further for more file types.
class open_gfa():
    def __init__(self, path, mode):
        self.path = path
        self.mode = mode

    def __enter__(self):
        self.fp = self._gfa_helper(self.path, self.mode)
        return self.fp

    def __exit__(self, type, value, traceback):
        self.fp.close()

    @staticmethod
    def _gfa_helper(file_path, mode):
        file_type = eval_file_type(file_path)
        if file_type:
            if "gz" == file_type:
                return gzip.open(file_path, 'r')
            else:
                raise Exception
        else:
            return open(file_path, 'r')


@profile
def read_gfa(file_path, haplotype=False):
    """
    Read GFA file, parse and return graph
    """

    graph = nx.DiGraph()

    if haplotype:
        haplotypes = get_haplo_header(file_path)
        graph.htypes = set(haplotypes)
        graph.hdict = dict([(str(i), j) for i, j in enumerate(haplotypes)])
    else:
        graph.htypes = 0

    with open_gfa(file_path, 'r') as file_in:
        for line in file_in:
            if line.startswith('S'):
                add_node(graph, line)
            elif line.startswith('L'):
                create_connection(graph, line)

    return graph


# A temporary alternative for 1000G data
@profile
def read_gfa_alt(file_path, haplotype=False):
    """
    Read GFA file, parse and return graph
    """

    graph = nx.DiGraph()

    if haplotype:
        haplotypes = get_haplo_header(file_path)
        graph.htypes = set(haplotypes)
        graph.hdict = dict([(str(i), j) for i, j in enumerate(haplotypes)])
    else:
        graph.htypes = 0

    with open_gfa(file_path, 'r') as file_in:
        for line in file_in:
            if line.startswith('S'):
                add_node_alt(graph, line)
            elif line.startswith('L'):
                create_connection(graph, line)

    return graph


@profile
def read_gfa_edge(file_path):
    graph = nx.DiGraph()
    haplotypes = get_haplo_header(file_path)
    graph.htypes = set(haplotypes)
    graph.hdict = dict([(str(i), j) for i, j in enumerate(haplotypes)])

    with open_gfa(file_path, 'r') as file_in:
        for line in file_in:
            if line.startswith('S'):
                add_node_alt(graph, line)
            elif line.startswith('L'):
                create_connection_alt(graph, line)

    return graph


def read_fasta(path, interval):
    with open(path, 'r') as f:
        while True:
            # TODO: use iter(partial()) instead of islice?
            n_lines = list(islice(f, interval))
            if not n_lines:
                break
            n_lines = [i.rstrip() for i in n_lines]
            yield n_lines


def all_paths(graph):
    """
    Return an iterator for all paths in a graph (concatenation of all sequences on nodes)
    Assumes a fully connected graph (i.e each node is reachable from any other node)
    """

    top_list = nx.topological_sort(graph)
    if top_list:
        for path in nx.all_simple_paths(graph, source=top_list[0], target=top_list[-1]):
            yield ''.join([graph.node[node]['sequence'] for node in path])


def get_kmers(seq, k, k_dict=None):
    """
    Return dictionary of k-mers for a given sequence and k, which may be an empty dict or one that is being filled
    """

    if k_dict is None:
        k_dict = defaultdict(int)

    for i in xrange(len(seq) - k + 1):
        k_dict[seq[i:i + k]] += 1
    return k_dict


def count_values(in_list):
    """
    Return the count of the number of elements in an arbitrarily nested list
    """

    # Deep copy of in_list
    if type(in_list) is dict:
        mod_list = list(in_list.values())
    else:
        mod_list = list(in_list)

    count = 0
    while mod_list:
        entry = mod_list.pop()
        if isinstance(entry, list):
            mod_list += entry
        else:
            count += 1
    return count


def concat_sequence(graph, node):
    """
    Return the concatenation of all sequence on a node, i.e:
        return (prefix + sequence + suffix)
    """

    seq_attributes = sorted(
        ATTRIBUTE_SET.intersection(graph.node[node].keys()))

    return ''.join([graph.node[node][i] for i in seq_attributes])


def get_stats(graph, k):
    """
    Number of nodes
    Number of edges
    Number of nodes with prefix
    Number of nodes with suffix
    Number of nodes with both a suffix and prefix
    Number of nodes without a prefix and suffix
    Average sequence length per node
    Shortest sequence length over all nodes
    Longest sequence length over all nodes
    Number of splits
    Numbers of collapses
    Haplotypes
    Number of haplotypes
    """

    num_nodes = graph.number_of_nodes()

    nodes = graph.nodes()
    nodes_prefix = [i for i in nodes if 'prefix' in graph.node[i].keys()]
    nodes_suffix = [i for i in nodes if 'suffix' in graph.node[i].keys()]

    nodes_both = set(nodes_prefix).intersection(nodes_suffix)
    nodes_prefix = set(nodes_prefix).difference(nodes_both)
    nodes_suffix = set(nodes_suffix).difference(nodes_both)
    nodes_none = set(nodes).difference(nodes_both, nodes_suffix, nodes_prefix)

    avg_len = len(concat_sequence(graph, nodes[0]))
    min_len, max_len = avg_len, avg_len
    for node in nodes[1:]:
        cur_len = len(concat_sequence(graph, node))
        avg_len += cur_len
        max_len = max(max_len, cur_len)
        min_len = min(min_len, cur_len)
    avg_len /= float(num_nodes)

    stat_str = "Graph statistics for k = %i:\n\
    \tNumber of nodes: %i\n\
    \tNumber of edges: %i\n\
    \tNumber of nodes with prefix: %i, suffix: %i, both: %i, or neither: %i\n\
    \tMinimum sequence length: %i\n\
    \tMaximum sequence length: %i\n\
    \tAverage sequence length: %.2f\
    " % (k, num_nodes, graph.number_of_edges(), len(nodes_prefix), len(nodes_suffix), len(nodes_both), len(nodes_none), min_len, max_len, avg_len)

    if hasattr(graph, 'n_split') and graph.n_split > 0:
        stat_str += "\n\tNumber of splits: %i" % graph.n_split

    if hasattr(graph, 'n_collapse') and graph.n_collapse > 0:
        stat_str += "\n\tNumber of collapses: %i" % graph.n_collapse

    if graph.htypes:
        stat_str += "\n\tNumber of haplotypes: %i\n\
    \tHaplotypes: " % len(graph.htypes)
        stat_str += ' '.join(sorted(graph.htypes)[:5])

    return stat_str


def over_paths(graph, k, node):
    """
    Helper function of get_k_paths()
    Yields all k-length paths over nodes given a starting node
    """

    extend_queue = deque()
    map(extend_queue.append, [
        [(succ, len(graph.node[succ]['sequence']))] for succ in graph.successors(node)])

    while extend_queue:
        extend = extend_queue.popleft()
        last_node, last_seq_len = extend[-1]
        if k - 1 > last_seq_len:
            for succ in graph.successors(last_node):
                tmp = extend[:]
                tmp.append(
                    (succ, last_seq_len + len(graph.node[succ]['sequence'])))
                extend_queue.append(tmp)

        yield "".join([graph.node[ex_node]['sequence'] for (ex_node, seq_len) in extend])[:k - 1]


def get_k_paths(graph, k, k_dict=None, source=None, sink=None):
    """
    Return all kmers in a graph between a (source, sink) both within nodes and across nodes.
    """

    if k_dict is None:
        k_dict = defaultdict(int)

    if source is None:
        node_queue = deque(
            [node for node, degree in graph.in_degree().items() if degree == 0])
    else:
        node_queue = deque([source])

    seen_set = set()

    if sink is not None:
        # All nodes that come after the sink should not be evaluated
        for s, t in nx.dfs_edges(graph, sink):
            seen_set.add(t)

    # Continue k-mer extraction until queue is empty
    while node_queue:
        node = node_queue.popleft()

        node_succ = set(graph.successors(node))

        # This is the costliest operation since seen_set will keep growing...
        unseen_set = node_succ - seen_set

        # Only add successors to the queue that are not in seen_set
        map(node_queue.append, unseen_set)
        seen_set.update(unseen_set)

        seq = graph.node[node]['sequence']

        # k_min denotes the index range we can extract k-mers from [0, k_min]
        # It also gives the starting point where we need to find k-mers across nodes up until k_max
        # In the event that k_min is < 0, k does not fit in the sequence and we start at index 0
        # (k_max is unused for now)
        # (k_min, k_max) = max(0, len(seq) - k + 1), len(seq)
        k_min = max(0, len(seq) - k + 1)

        # K-mers that fall within the current node
        for i in xrange(k_min):
            kmer = seq[i:i + k]
            k_dict[kmer] += 1

        # K-mers across paths
        for path in over_paths(graph, k, node):
            n_seq = seq[k_min:] + path
            get_kmers(n_seq, k, k_dict)

    return k_dict


def consolidate_ranges(r1, r2):
    return map(lambda x: x + (max(r1) + 1), r2)


def get_haplotype_k_paths(graph, k, k_dict=None, n_bits=2):
    """
    Extract all k-length haplotype supported paths in a DAG returning a <kmer,count> dictionary
    """

    if k_dict is None:
        k_dict = defaultdict(int)

    bit_length = len(graph.htypes) * n_bits
    haplotype_sets = [x() for x in [lambda m=i * n_bits: set([m, m + 1])
                                    for i, j in enumerate(sorted(graph.htypes))]]

    # K-mers contained within nodes
    for node in graph.nodes():
        get_kmers(graph.node[node]['sequence'], k, k_dict)

    # K-mers across nodes
    for haplotype in haplotype_sets:
        haplotype = bitarray(
            ''.join(['1' if i in haplotype else '0' for i in xrange(bit_length)]))
        for node in graph.nodes():
            edge_list = [(u, v) for u, v in graph.out_edges(
                node) if any(haplotype & graph[u][v])]

            if edge_list:
                # Build the base sequence from the sequence on the current node. this is the part where k-mer
                # extraction requires extension across nodes, e.g for k = 5
                # Node sequence: ATTAGAGATAGAGAGGGAAG
                #                                ****
                # Base sequence:                 GAAG
                add_seq = [graph.node[node]['sequence'][
                    max(0, len(graph.node[node]['sequence']) - k + 1):]]

                # The remainder is the maximum number of bases that should be added to extract a k-mer starting
                # from the last position of the base sequence (0th element in
                # add_seq)
                remainder = k - 1

                while (remainder > 0 and edge_list):
                    (_, v), = edge_list
                    len_seq = len(graph.node[v]['sequence'])
                    add_seq.append(graph.node[v][
                                   'sequence'][:remainder])
                    remainder -= len_seq

                    # As long as there is a remainder try to get the next edge
                    if remainder:
                        edge_list = [(u_out, v_out) for u_out, v_out in graph.out_edges(
                            v) if any(haplotype & graph[u_out][v_out])]

                get_kmers(''.join(add_seq), k, k_dict)

    return k_dict


def move_haplotype_to_edges(graph):
    """
    Move haplotypes encoded on nodes to the edges
    """

    for node in nx.topological_sort(graph):
        for pred in graph.predecessors(node):
            u_haplo = graph.node[pred]['haplotype']
            v_haplo = graph.node[node]['haplotype']

            graph[pred][node] = u_haplo & v_haplo
            graph.node[pred]['haplotype'] = ~(~u_haplo | v_haplo)

    for node in graph.nodes():
        del graph.node[node]['haplotype']


@profile
def set_interval_nodes(graph):
    """
    Initialize intervals on all nodes of the graph
    """

    for node in graph.nodes():
        graph.node[node]['sequence_from'] = [NodeInterval(
            node, (0, len(graph.node[node]['sequence'])))]


class NodeInterval(object):
    def __init__(self, node, interval):
        self.node = node
        self.interval = interval

    @property
    def node(self):
        return self._node

    @node.setter
    def node(self, n):
        if n is None:
            raise Exception("Node cannot be empty")
        self._node = n

    @property
    def interval(self):
        return self._interval

    @interval.setter
    def interval(self, i):
        if not isinstance(i, tuple):
            raise Exception("Interval must be a tuple")
        self._interval = i


# TODO: generalize
def interval_ins_tail(node_from, node_to, new=False):
    # node_from = deepcopy(node_from)  # TODO: Is this still necessary?

    if new:
        node_to_iter = deepcopy(node_to)
    else:
        node_to_iter = node_to

    for i in node_from:
        node_to_iter.append(i)

    if new:
        return node_to_iter


# TODO: generalize
def interval_ins_head(node_from, node_to, new=False):
    # node_from = deepcopy(node_from)  # TODO: Is this still necessary?

    if new:
        node_to_iter = deepcopy(node_to)
    else:
        node_to_iter = node_to

    node_to_iter[:0] = node_from

    if new:
        return node_to_iter


def simple_suffix(node_from, k_1):
    remainder = k_1
    min_interval = 0
    max_interval = 0
    suffix_interval = []

    for i in node_from:
        interval_len = i.interval[1]
        remainder -= interval_len

        # The entire interval fits
        if remainder > 0:
            max_interval += interval_len
            suffix_interval.append(NodeInterval(
                i.node, (0, max_interval - min_interval)))
            min_interval = max_interval
            max_interval = min_interval
        # Part of the current interval is required
        else:
            max_interval = remainder + interval_len
            suffix_interval.append(NodeInterval(
                i.node, (0, max_interval)))
            break

    return suffix_interval


def simple_prefix(node_from, k_1, graph=None):
    remainder = k_1
    min_interval = 0
    max_interval = k_1
    prefix_interval = []

    for i in reversed(node_from):
        interval_len = i.interval[1]
        remainder -= interval_len

        # The entire interval fits
        if remainder > 0:
            min_interval = max_interval - interval_len
            prefix_interval.insert(0, NodeInterval(
                i.node, (0, max_interval - min_interval)))
            max_interval = min_interval
        # Part of the current interval is required
        else:
            min_interval = interval_len - max_interval
            max_interval = interval_len
            prefix_interval.insert(0, NodeInterval(
                i.node, (min_interval, max_interval)))
            break

    return prefix_interval


def emit_identifier(interval):
    return '_'.join([str(i.node) + "+{0},{1}".format(*i.interval) for i in interval])


def generate_identifier(graph, node):
    node_attributes = graph.node[node]
    sequence_interval = node_attributes['sequence_from']
    if node_attributes.get('prefix_from'):
        sequence_interval = interval_ins_head(
            node_attributes['prefix_from'], sequence_interval, True)

    if node_attributes.get('suffix_from'):
        sequence_interval = interval_ins_tail(
            node_attributes['suffix_from'], sequence_interval, True)

    return '_'.join([str(i.node) + "+{0},{1}".format(*i.interval) for i in sequence_interval])


def write_fasta(graph, path=None):
    if not graph.nodes():
        exit(1)

    if path:
        with open(path, 'w') as f:
            logger.debug("Start writing to {}".format(path))
            for node in graph.nodes():
                ident = generate_identifier(graph.node[node])
                logger.debug("Writing node {}".format(node))
                f.write("{}\n{}\n".format(ident, concat_sequence(graph, node)))
        logger.debug("Finished writing to {}".format(path))
    else:
        for node in graph.nodes():
            print "{}\t{}\n{}".format(node, generate_identifier(graph.node[node]), concat_sequence(graph, node))


def validate_interval_graph(graph, a_graph):
    for node in graph.nodes():
        identifier = generate_identifier(graph.node[node])
        node_sequence = concat_sequence(graph, node)

        intervals = [i.split(',') for i in identifier.split(';')]
        for i in intervals:
            i_node = int(i[0].lstrip('>'))
            i_li = int(i[1].lstrip('<'))
            i_ri = int(i[2].rstrip('>'))
            i = (i_node, i_li, i_ri)
            if not node_sequence[i[1]:i[2]] in a_graph.node[i[0]]['sequence']:
                return False
        else:
            return True


def unw_longest_path(in_graph, destructive=False):
    """
    Input DAG w/ sequences on all nodes
    Return length of longest path
    """

    if destructive:
        graph = in_graph
    else:
        graph = deepcopy(in_graph)

    top_sort = nx.topological_sort(graph)

    graph.node[top_sort[0]]['dist'] = len(graph.node[top_sort[0]]['sequence'])

    for node in top_sort[1:]:
        node_predecessors = graph.predecessors(node)
        if len(node_predecessors) == 1:
            (pred,) = node_predecessors
            graph.node[node]['dist'] = graph.node[pred][
                'dist'] + len(graph.node[node]['sequence'])
        else:
            node_dist_list = [(graph.node[pred][
                               'dist'] + len(graph.node[node]['sequence'])) for pred in node_predecessors]
            graph.node[node]['dist'] = max(node_dist_list)

    return graph.node[top_sort[-1]]['dist']


def get_complement(seq):
    complement = maketrans('ATCG', 'TAGC')
    return seq.translate(complement)


def graph_to_dot(graph):
    for node in nx.topological_sort(graph):
        print '{} {}"];'.format(node, '[label="' + graph.node[node]['sequence'])
        successors = graph.successors(node)
        if successors:
            for succ in successors:
                print "{};".format(' -> '.join(map(str, [node, succ])))


def get_total_sequence_fq(fq_path, interval=4, read_length=None):
    # Length of reads are variable, required to iterate over all sequence and
    # determine length
    if read_length is None:
        line_count = 0
        sequence_count = 0
        with open(fq_path, 'r') as f:
            interval_line = islice(f, 1, None, interval)
            for line in interval_line:
                line_count += 1
                sequence_count += len(line.strip())
        return ((sequence_count / float(line_count)), sequence_count)

    # Length of reads is unknown but it is the same for all reads, do a peek at the 2nd line afterwards do the same as
    # when read_length is known
    elif read_length is True:
        with open(fq_path, 'r') as f:
            for i, line in enumerate(f):
                if i == 1:
                    read_length = len(line.strip())
                    break
        return (read_length, (get_n_of_lines(fq_path) / interval) * read_length)

    # Length of reads is known beforehand
    elif isinstance(read_length, int):
        return (read_length, (get_n_of_lines(fq_path) / interval) * read_length)


def get_haplotype_paths(graph):
    haplotype_dict = defaultdict(str)

    for node in nx.topological_sort(graph):
        for haplotype in graph.node[node]['haplotype']:
            haplotype_dict[haplotype] = ''.join(
                [haplotype_dict[haplotype], graph.node[node]['sequence']])
    return haplotype_dict


def make_sequence(graph, node, prefix=''):
    if prefix:
        return ">{0}|{1}\n{2}\n".format(prefix, generate_identifier(graph, node), concat_sequence(graph, node))
    else:
        return ">{0}\n{1}\n".format(generate_identifier(graph, node), concat_sequence(graph, node))


def go(node, graph, radius, direction_edges, index):
    if radius > 0:
        edges_to_mark = direction_edges(node)
        for neighbor_edge in edges_to_mark:
            for i in neighbor_edge:
                graph.node[i]['m'] = True
            go(neighbor_edge[index], graph, radius - 1, direction_edges, index)


def mark_radius(graph, edge, radius=1):
    # Merge these two
    map(partial(go, graph=graph, radius=radius,
                direction_edges=graph.in_edges, index=0), edge)
    map(partial(go, graph=graph, radius=radius,
                direction_edges=graph.out_edges, index=1), edge)


def is_marked_edge(graph, edge):
    for node in edge:
        node_data = graph.node.get(node, True)
        # TODO: Why not just this?
        # return isinstance(node_data, bool) or node_data.get('m', True)
        if isinstance(node_data, bool):
            return True
        if node_data.get('m', True):
            return True
    return False


def connected_component_subgraphs(graph, s_buffer, prefix=''):
    for component in get_connected_components_edit(graph):
        if len(component) > 1:
            yield graph.subgraph(component).copy()
            graph.remove_nodes_from(component)
        else:
            node = component.pop()
            s_buffer.append(make_sequence(graph, node, prefix))
            graph.remove_node(node)


def bfs_directed(graph, source):
    seen_set = set()
    current_step = set([source])

    while current_step:
        next_step = current_step
        current_step = set()
        for node in next_step:
            if node not in seen_set:
                yield node
                seen_set.add(node)
                current_step.update(graph.succ[node])
                current_step.update(graph.pred[node])


def get_connected_components(graph):
    seen_set = set([])
    for node in graph:
        if node not in seen_set:
            component = set(bfs_directed(graph, node))
            yield component
            seen_set.update(component)


def get_connected_components_edit(graph):
    """
    Alternate variant of get_connected_components that allows the graph to be modified
    """

    seen_set = set([])
    graph_nodes = list(graph.nodes())
    for node in graph_nodes:
        if node not in seen_set:
            component = set(bfs_directed(graph, node))
            yield component
            seen_set.update(component)


def write_gfa(graph, path, samples):
    samples = list(samples)
    sample_lookup = dict((j, str(i)) for i, j in enumerate(samples))
    header = "H\tVN:Z:1.0\nH\tORI:{0}\n".format(';'.join(samples))
    with open(path, 'wb') as out_f:
        out_f.write(header)
        for idx, node in enumerate(graph.nodes(), start=1):
            if idx % 10000 == 0:
                logger.info("Processed {} nodes".format(idx + 1))

            # Segment (node)
            out_f.write("S\t{0}\t{1}\t*\tORI:Z:{2}\n".format(node, re.sub(
                'N+', 'N', graph.node[node]['sequence']), ';'.join([sample_lookup[i] for i in graph.node[node]['haplotype']])))

            # Link (edge)
            for successor in graph.successors(node):
                out_f.write("L\t{0}\t+\t{1}\t+\t0M\n".format(node, successor))


def is_valid_interval(path, mod_graph=None, input_graph=None, node=None):
    """
    Return path validity
    """

    if mod_graph:
        s = []

    start, end = path[0].interval
    if start < 0 or start > end:
        return False

    cur_node = path[0].node

    if mod_graph:
        s.append(input_graph.node[cur_node]['sequence'][start:end])

    for entry in path[1:]:
        start, end = entry.interval

        if start != 0 or end <= 0 or cur_node >= entry.node:
            return False
        cur_node = entry.node
        if mod_graph:
            s.append(input_graph.node[cur_node]['sequence'][start:end])

    if mod_graph:
        return ''.join(s) == mod_graph.node[node]['sequence']

    return True


######################################################################
# ooo.                                              o              8 #
# 8  `8.                                            8              8 #
# 8   `8 .oPYo. .oPYo. oPYo. .oPYo. .oPYo. .oPYo.  o8P .oPYo. .oPYo8 #
# 8    8 8oooo8 8    8 8  `' 8oooo8 8    ' .oooo8   8  8oooo8 8    8 #
# 8   .P 8.     8    8 8     8.     8    . 8    8   8  8.     8    8 #
# 8ooo'  `Yooo' 8YooP' 8     `Yooo' `YooP' `YooP8   8  `Yooo' `YooP' #
#               8                                                    #
#               8                                                    #
######################################################################


def stretch_gaps(graph):
    """
    Split nodes with gap(s) (N) into new nodes that either contain gap sequence or 'proper' sequence
    """

    graph.max_id = max([i for i in graph.nodes()])

    # Iterate over all nodes that have at least one gap in the sequence
    for node in graph.nodes():
        if 'N' in graph.node[node]['sequence']:
            # Split the sequence in a list based on gaps, retain the delimiter
            seq_to_split = filter(None, re.split(
                '([N]*)', graph.node[node]['sequence']))
            seq_to_split_len = len(seq_to_split)

            if seq_to_split_len > 1:
                orig_out_edges = graph.out_edges(node)
                orig_successors = graph.successors(node)

                # Process first node
                graph.remove_edges_from(orig_out_edges)
                graph.node[node]['sequence'] = seq_to_split[0]

                # Process nodes between first and last nodes
                prev_node = node
                for i in seq_to_split[1:-1]:
                    graph.max_id += 1
                    graph.add_node(graph.max_id)
                    graph.node[graph.max_id]['origin'] = node
                    graph.node[graph.max_id]['sequence'] = i
                    graph.add_edge(prev_node, graph.max_id)

                    prev_node = graph.max_id

                # Process last node
                graph.max_id += 1
                graph.add_node(graph.max_id)
                graph.node[graph.max_id]['origin'] = node
                graph.node[graph.max_id]['sequence'] = seq_to_split[-1]
                graph.add_edge(prev_node, graph.max_id)
                graph.add_edges_from([(graph.max_id, i)
                                      for i in orig_successors])


def stretch_gaps_cut_excess(graph, node, seq_to_split):
    orig_out_edges = graph.out_edges(node)
    orig_in_edges = graph.in_edges(node)
    orig_successors = graph.successors(node)

    # Process first node
    graph.remove_edges_from(orig_out_edges)
    graph.node[node]['sequence'] = seq_to_split[0]
    if 'N' in graph.node[node]['sequence']:
        graph.remove_edges_from(orig_in_edges)

    # Process nodes between first and last
    prev_node = node
    for i in seq_to_split[1:-1]:
        graph.max_id += 1
        graph.add_node(graph.max_id)
        graph.node[graph.max_id]['origin'] = node
        graph.node[graph.max_id]['sequence'] = i

        if 'haplotype' in graph.node[node]:
            graph.node[graph.max_id][
                'haplotype'] = graph.node[node]['haplotype']

        if 'N' not in i and 'N' not in graph.node[prev_node]['sequence']:
            graph.add_edge(prev_node, graph.max_id)

        prev_node = graph.max_id

    # Process last node
    graph.max_id += 1
    graph.add_node(graph.max_id)
    graph.node[graph.max_id]['origin'] = node
    graph.node[graph.max_id]['sequence'] = seq_to_split[-1]

    if 'haplotype' in graph.node[node]:
        graph.node[graph.max_id]['haplotype'] = graph.node[node]['haplotype']

    if 'N' not in seq_to_split[-1]:
        graph.add_edges_from([(graph.max_id, i) for i in orig_successors])
        if 'N' not in graph.node[prev_node]['sequence']:
            graph.add_edge(prev_node, graph.max_id)


def stretch_gaps_cut(graph):
    """
    Split nodes with gap(s) (N) into new nodes that either contain gap sequence or 'proper' sequence
    Nodes that are adjacent to gaps should not be extended to them OR over them, so we can safely remove them!
    """

    graph.max_id = max([i for i in graph.nodes()])

    # Iterate over all nodes and find those that have at least one gap in the
    # sequence
    for node in graph.nodes():
        if 'N' in graph.node[node]['sequence']:
            # Split the sequence in a list based on gaps, retain the delimiter
            # and remove any None sequences
            seq_to_split = filter(None, re.split(
                '([N]*)', graph.node[node]['sequence']))
            seq_to_split_len = len(seq_to_split)

            # Easy case: there is one entry, which must be an N. Just cut the
            # in/out-going edges
            if seq_to_split_len == 1:
                graph.remove_node(node)
            # There is more sequence meaning we have to split things up
            else:
                stretch_gaps_cut_excess(graph, node, seq_to_split)

    for node in graph.nodes():
        if 'N' in graph.node[node]['sequence']:
            graph.remove_node(node)


def collapse_chains(graph):
    """
    Test whether collapsing chains after the indexing 'breaks' will still mean we capture all k-paths
    """

    starts_of_chain = [node for node, pred in graph.pred.items() if len(
        pred) == 0 and len(graph.successors(node)) > 0]
    for node in starts_of_chain:
        successor_dict = nx.dfs_successors(graph, node)
        seq = graph.node[node]['sequence']
        while True:
            try:
                tmp = node
                graph.remove_node(node)
                (node,) = successor_dict[tmp]
                seq = ''.join([seq, graph.node[node]['sequence']])
            except:
                break
        graph.max_id += 1
        graph.add_node(graph.max_id, sequence=seq)


if __name__ == '__main__':

    print len(get_kmers("ATCCTAGCTTCCGCCCAC", 4).keys())
    # graph = read_gfa("/Users/tomb/Downloads/nvc-filt.gfa")

    # topo = nx.topological_sort(graph)

    # sPaths = 1

    # for node in topo:
    #     sPaths *= max(1, graph.out_degree(node))

    # print sPaths

    # graph = nx.DiGraph()
    # graph.htypes = set(['1', '2', '3'])
    # #
    # graph.add_node(0, sequence='TACG', haplotype=set(['1', '2', '3']))
    # # | 1
    # graph.add_node(1, sequence='G', haplotype=set(['1', '3']))
    # graph.add_node(2, sequence='C', haplotype=set(['2']))
    # # | 2
    # graph.add_node(3, sequence='C', haplotype=set(['1', '2', '3']))
    # # | 3
    # graph.add_node(4, sequence='C', haplotype=set(['1']))
    # graph.add_node(5, sequence='T', haplotype=set(['3']))
    # graph.add_node(6, sequence='GTT', haplotype=set(['2']))
    # # | 4
    # graph.add_node(7, sequence='GG', haplotype=set(['1', '2', '3']))

    # #
    # graph.add_edge(0, 1)
    # graph.add_edge(0, 2)
    # #
    # graph.add_edge(1, 3)
    # graph.add_edge(2, 3)
    # #
    # graph.add_edge(3, 4)
    # graph.add_edge(3, 5)
    # graph.add_edge(3, 6)
    # #
    # graph.add_edge(4, 7)
    # graph.add_edge(5, 7)
    # graph.add_edge(6, 7)

    # write_gfa(graph, "/Users/tomb/Desktop/pipeline.gfa", graph.htypes)

    # print get_kmers('TACGGCGCTGG', 3).keys()
