from collections import deque
from copy import deepcopy
from ctypes import c_char_p
from graph_functions import profile
from random import choice, getrandbits
from time import sleep
from uuid import uuid4
import argparse
import colorer
import graph_functions as gf
import logging
import multiprocessing
import os
import signal
import sys

logger = logging.getLogger(__name__)


try:
    import cPickle as pickle
except ImportError:
    import pickle


def signal_kill(signal, frame):
    print 'Killed!'
    sys.exit(0)


signal.signal(signal.SIGINT, signal_kill)
logger = logging.getLogger(__name__)

MIN_GRAPH_SIZE = 500
CC_CUTOFF_SIZE = 250
MARK_RADIUS = 1
ATTRIBUTE_SET = set(['prefix', 'suffix', 'prefix_from', 'suffix_from'])


def subsequence_extension(graph, edge, simple_node_idx):
    """
    Attempt to extend subsequence between a pair of nodes
    """

    cur_node = edge[simple_node_idx]
    affix_node = edge[simple_node_idx ^ 1]

    affix_node_sequence = graph.node[affix_node]['sequence']

    if simple_node_idx:
        direction_label = 'prefix'
    else:
        direction_label = 'suffix'

    # Sufficient sequence, safe to extend
    if len(affix_node_sequence) >= graph.k_extend:
        # Prefix
        if simple_node_idx:
            logger.debug("EXTENSION: {} got prefix from {}".format(
                cur_node, affix_node))

            graph.node[cur_node][
                'prefix'] = affix_node_sequence[-graph.k_extend:]
            graph.node[cur_node]['prefix_from'] = gf.simple_prefix(
                graph.node[affix_node]['sequence_from'], graph.k_extend, graph)

            logger.debug("  * Prefix {} build from simple_prefix({}, {}, graph)".format(gf.emit_identifier(graph.node[
                cur_node]['prefix_from']), gf.emit_identifier(graph.node[affix_node]['sequence_from']), graph.k_extend))
            if logger.getEffectiveLevel() == 10 and not gf.is_valid_interval(graph.node[cur_node]['prefix_from']):
                logger.debug("  * Prefix: invalid - {} build from simple_prefix({}, {}, graph)".format(gf.emit_identifier(graph.node[
                    cur_node]['prefix_from']), gf.emit_identifier(graph.node[affix_node]['sequence_from']), graph.k_extend))
                exit(1)

        # Suffix
        else:
            logger.debug("EXTENSION: {} got suffix from {}".format(
                cur_node, affix_node))

            graph.node[cur_node][
                'suffix'] = affix_node_sequence[:graph.k_extend]
            graph.node[cur_node]['suffix_from'] = gf.simple_suffix(
                graph.node[affix_node]['sequence_from'], graph.k_extend)

            logger.debug("  * Suffix {} build from simple_suffix({}, {}, graph)".format(gf.emit_identifier(graph.node[
                cur_node]['suffix_from']), gf.emit_identifier(graph.node[affix_node]['sequence_from']), graph.k_extend))
            if logger.getEffectiveLevel() == 10 and not gf.is_valid_interval(graph.node[cur_node]['suffix_from']):
                logger.debug("  * Suffix: invalid - {} build from simple_suffix({}, {}, graph)".format(gf.emit_identifier(graph.node[
                    cur_node]['suffix_from']), gf.emit_identifier(graph.node[affix_node]['sequence_from']), graph.k_extend))
                exit(1)

        logger.debug("  * Deleting edge {}".format(edge))
        graph.remove_edge(*edge)
        return True
    else:
        # If the other node has been extended we can extend across it
        if direction_label in graph.node[affix_node]:
            if simple_node_idx:
                logger.debug("EXTENSION: {} got prefix from and across {}".format(
                    cur_node, affix_node))

                graph.node[cur_node]['prefix'] = graph.node[affix_node][
                    'prefix'][-(graph.k_extend - len(affix_node_sequence)):] + affix_node_sequence
                in_node = graph.node[affix_node]['sequence_from']
                over_node = gf.simple_prefix(graph.node[affix_node][
                                             'prefix_from'], graph.k_extend - len(affix_node_sequence), graph)

                logger.debug("  * Prefix across: {} build from simple_prefix({}, {})".format(gf.emit_identifier(over_node),
                                                                                             gf.emit_identifier(graph.node[affix_node]['prefix_from']), graph.k_extend - len(affix_node_sequence)))

                graph.node[cur_node]['prefix_from'] = gf.interval_ins_tail(
                    in_node, over_node, True)

                logger.debug("  * Prefix across {} build from interval_ins_tail({}, {}, True)".format(gf.emit_identifier(
                    graph.node[cur_node]['prefix_from']), gf.emit_identifier(in_node), gf.emit_identifier(over_node)))
                if logger.getEffectiveLevel() == 10 and not gf.is_valid_interval(graph.node[cur_node]['prefix_from']):
                    logger.debug("  * Prefix across: invalid - {} build from interval_ins_tail({}, {}, True)".format(gf.emit_identifier(
                        graph.node[cur_node]['prefix_from']), gf.emit_identifier(in_node), gf.emit_identifier(over_node)))
                    exit(1)

            else:
                logger.debug("EXTENSION: {} got suffix from and across {}".format(
                    cur_node, affix_node))

                graph.node[cur_node]['suffix'] = affix_node_sequence + graph.node[
                    affix_node]['suffix'][:graph.k_extend - len(affix_node_sequence)]
                in_node = graph.node[affix_node]['sequence_from']
                over_node = gf.simple_suffix(graph.node[affix_node][
                    'suffix_from'], graph.k_extend - len(affix_node_sequence))

                logger.debug("  * Suffix across; Before simple_suffix: {}; After simple_suffix: {}".format(
                    gf.emit_identifier(graph.node[affix_node]['suffix_from']), gf.emit_identifier(over_node)))

                graph.node[cur_node]['suffix_from'] = gf.interval_ins_head(
                    in_node, over_node, True)

                logger.debug("  * Suffix across: {} build from interval_ins_head({}, {}, True)".format(gf.emit_identifier(
                    graph.node[cur_node]['suffix_from']), gf.emit_identifier(in_node), gf.emit_identifier(over_node)))
                if logger.getEffectiveLevel() == 10 and not gf.is_valid_interval(graph.node[cur_node]['suffix_from']):
                    logger.debug("  * Suffix across: invalid - {} build from interval_ins_head({}, {}, True)".format(gf.emit_identifier(
                        graph.node[cur_node]['suffix_from']), gf.emit_identifier(in_node), gf.emit_identifier(over_node)))
                    exit(1)

            graph.remove_edge(*edge)
            logger.debug("  * Removing edge {}".format(edge))
            return True
        else:
            return False


def node_collapser(graph, node, edge_idx):
    """
    Attempt to collapse a pair of nodes into a single node
    """

    if edge_idx:
        (neighbor_node,) = graph.predecessors(node)
    else:
        (neighbor_node,) = graph.successors(node)

    assert(len(graph.successors(neighbor_node)
               if edge_idx else graph.predecessors(neighbor_node)) == 1)

    if edge_idx:
        edge = (neighbor_node, node)
    else:
        edge = (node, neighbor_node)

    for attribute in [i for i in graph.node[neighbor_node].keys() if i in ATTRIBUTE_SET]:
        graph.node[node][attribute] = deepcopy(
            graph.node[neighbor_node][attribute])

    if edge_idx:
        logger.debug("COLLAPSE: {} --> {}".format(
            neighbor_node, node))

        graph.node[node]['sequence'] = graph.node[neighbor_node][
            'sequence'] + graph.node[node]['sequence']
        inherit_edges_list = graph.in_edges(neighbor_node)
        gf.interval_ins_head(graph.node[neighbor_node][
            'sequence_from'], graph.node[node]['sequence_from'])
    else:
        logger.debug("COLLAPSE: {} --> {}".format(
            neighbor_node, node))

        graph.node[node][
            'sequence'] += graph.node[neighbor_node]['sequence']
        inherit_edges_list = graph.out_edges(neighbor_node)
        gf.interval_ins_tail(graph.node[neighbor_node][
            'sequence_from'], graph.node[node]['sequence_from'])

    for edge in inherit_edges_list:
        if edge_idx:
            new_edge = (edge[edge_idx ^ 1], node)
        else:
            new_edge = (node, edge[edge_idx ^ 1])

        graph.add_edge(*new_edge)

        if graph.htypes:
            graph[new_edge[0]][new_edge[1]
                               ] = graph[edge[0]][edge[1]]
            logger.debug("  * Added edge {} with haplotype {}".format(
                new_edge, graph[new_edge[0]][new_edge[1]].to01()))

    graph.remove_node(neighbor_node)
    logger.debug("  * Delete node {}".format(neighbor_node))


def node_splitter(graph, node, edge_idx):
    """
    Split nodes to resolve edge ambiguities
    """

    if edge_idx:
        edge_list = graph.in_edges(node)
        inherit_edge_list = graph.out_edges(node)
    else:
        edge_list = graph.out_edges(node)
        inherit_edge_list = graph.in_edges(node)

    processed_nodes = [node]
    logger.debug("DUPLICATE: node {}".format(node))

    # Iterate over all but one of the in/out-coming edges of the original
    # node, for each introducing a duplicate
    for edge in edge_list[1:]:
        graph.max_id += 1
        graph.add_node(graph.max_id)
        logger.debug("  * Added node {}".format(graph.max_id))
        processed_nodes.append(graph.max_id)

        # Inherit the attributes of the original node
        for key in graph.node[node]:
            graph.node[graph.max_id][key] = deepcopy(
                graph.node[node][key])

        # Single edge taken from the original node
        if edge_idx:
            new_edge = (edge[edge_idx ^ 1], graph.max_id)
        else:
            new_edge = (graph.max_id, edge[edge_idx ^ 1])

        graph.add_edge(*new_edge)
        logger.debug("  * Added edge {}".format(
            new_edge))

        if graph.htypes:
            graph[new_edge[0]][new_edge[1]] = graph[edge[0]][edge[1]]
            logger.debug("  * Added edge {} with haplotype {}".format(
                new_edge, graph[new_edge[0]][new_edge[1]].to01()))
            # Get edges in the opposite direction that intersect with the
            # haplotypes of the newly duplicated node
            valid_edges = [i for i in inherit_edge_list if any(
                graph[i[0]][i[1]] & graph[new_edge[0]][new_edge[1]])]

            for i in valid_edges:
                # Introduce edges in the opposite direction
                if edge_idx:
                    new_opposite_edge = (graph.max_id, i[edge_idx])
                else:
                    new_opposite_edge = (i[edge_idx], graph.max_id)
                graph.add_edge(*new_opposite_edge)

                # In- and out-going haplotypes should be consistent
                graph[new_opposite_edge[0]][new_opposite_edge[1]] = graph[
                    i[0]][i[1]] & graph[new_edge[0]][new_edge[1]]
                logger.debug("  * Added edge {} with haplotype {}".format(
                    new_opposite_edge, graph[new_opposite_edge[0]][new_opposite_edge[1]].to01()))
        else:
            for i in inherit_edge_list:
                new_opposite_edge = (graph.max_id, i[edge_idx]) if edge_idx else (
                    i[edge_idx], graph.max_id)
                graph.add_edge(*new_opposite_edge)
                logger.debug("  * Added edge {}".format(new_opposite_edge))

    logger.debug("  * Deleting edges from: {}".format(edge_list[1:]))
    graph.remove_edges_from(edge_list[1:])

    if graph.htypes:
        for i in inherit_edge_list:
            # Remove edge if there is no co-occurence of haplotypes
            if not any(graph[edge_list[0][0]][edge_list[0][1]] & graph[i[0]][i[1]]):
                graph.remove_edge(*i)
                logger.debug("  * Delete edge: {}".format(i))
            # Ensure only haplotypes remain that actually enter the node
            else:
                graph[i[0]][i[1]] = graph[i[0]][i[1]] & graph[
                    edge_list[0][0]][edge_list[0][1]]
                logger.debug("  * Modified haplotype edge: {} to {}".format(
                    i, graph[i[0]][i[1]].to01()))


def get_sni(graph, edge, edge_degree):
    simple_node_idx = [idx for idx, val in enumerate(edge_degree) if val == 1]

    if len(simple_node_idx) == 1:
        (simple_node_idx,) = simple_node_idx

    # Resolve cases where both nodes are possible candidates, get best choice
    elif len(simple_node_idx) > 1:
        # See which of the nodes has sufficient sequence or is already extended
        candidate_nodes = [idx for idx in simple_node_idx if len(graph.node[edge[idx]][
                                                                 'sequence']) >= graph.k_extend or graph.node[edge[idx]].get('suffix' if idx else 'prefix')]
        # Get rid of the node that is a candidate for extension
        if len(candidate_nodes) == 1:
            simple_node_idx = [
                i for i in simple_node_idx if i not in candidate_nodes]

        if simple_node_idx:
            # If both nodes meet the requirements randomly pick one
            if len(simple_node_idx) > 1:
                (simple_node_idx,) = [choice(simple_node_idx)]
            else:
                (simple_node_idx,) = simple_node_idx
    else:
        return -1

    return simple_node_idx


def min_n_components(cc, n):
    count = 0
    for c in (len(component) > CC_CUTOFF_SIZE for component in iter(cc)):
        if c:
            count += 1
            if count >= n:
                return True
    return False


@profile
def first_sweep(graph):
    g_edges = graph.edges()
    for edge in g_edges:
        edge_degree = (gf.out_deg(graph, edge[0]), gf.in_deg(graph, edge[1]))
        simple_node_idx = get_sni(graph, edge, edge_degree)

        if simple_node_idx >= 0 and not all([i == 1 for i in edge_degree]):
            subsequence_extension(graph, edge, simple_node_idx)


def process_subgraph(graph):
    g_edges = graph.edges()

    # Mark all nodes as unmarked
    for node in graph.nodes():
        graph.node[node]['m'] = False

    # Process all edges that are unmarked
    for edge in g_edges:
        if not gf.is_marked_edge(graph, edge):
            edge_degree = (gf.out_deg(
                graph, edge[0]), gf.in_deg(graph, edge[1]))
            simple_node_idx = get_sni(graph, edge, edge_degree)

            if all([i == 1 for i in edge_degree]):
                gf.mark_radius(graph, edge, radius=MARK_RADIUS)
                rand_choice = getrandbits(1)
                node_collapser(graph, edge[rand_choice], rand_choice)
            elif simple_node_idx < 0:
                gf.mark_radius(graph, edge, radius=MARK_RADIUS)
                split_idx = min(idx for (idx, val)
                                in enumerate(edge_degree) if val > 1)
                node_splitter(graph, edge[split_idx], split_idx)
            else:
                isExtended = subsequence_extension(
                    graph, edge, simple_node_idx)
                if not isExtended:
                    gf.mark_radius(graph, edge, radius=MARK_RADIUS)
                    split_idx = min(idx for (idx, val)
                                    in enumerate(edge_degree) if val > 1)
                    node_splitter(graph, edge[split_idx], split_idx)


def full_index(graph, k, graph_queue, s_buffer, haplotype, prefix='', queue_counter=None):
    graph.k_extend = k - 1  #
    graph.htypes = haplotype

    graph.max_id = max([i for i in graph.nodes()])  #

    # Break apart the subgraph into smaller subgraphs
    if graph.number_of_edges() > MIN_GRAPH_SIZE:
        while graph.number_of_edges():
            process_subgraph(graph)

            co = []
            for component in list(gf.get_connected_components(graph)):
                if len(component) == 1:
                    node = next(iter(component))
                    s_buffer.append(gf.make_sequence(graph, node, prefix))
                    graph.remove_node(node)
                else:
                    co.append(component)

            # If the number of components exceeds the cutoff, split them into subgraphs
            if min_n_components(co, 2):
                extend_list = [graph.subgraph(component) for component in co]

                graph.clear()

                if extend_list:
                    if isinstance(graph_queue, deque):
                        graph_queue.extend(extend_list)
                    else:
                        for subgraph in extend_list:
                            graph_queue.put(subgraph)
                            queue_counter.increment()
    else:
        # Continously pass over the graph until all edges are resolved
        while graph.number_of_edges():
            process_subgraph(graph)


class Counter(object):
    def __init__(self, init=0):
        self.val = multiprocessing.Value('i', init)

    def increment(self, n=1):
        with self.val.get_lock():
            self.val.value += n

    def decrement(self, n=1):
        with self.val.get_lock():
            self.val.value -= n

    @property
    def value(self):
        return self.val.value


def worker_main(queue, counter, graph_counter, queue_counter, k, file_list, working_directory, uniq_id, haplotype, prefix=''):
    worker_name = multiprocessing.current_process().name
    logger.debug("Started: {}".format(worker_name))

    prefix_sequence = prefix.value

    out_path = os.path.join(
        working_directory, "{0}_{1}.tmp".format(worker_name, uniq_id))

    file_list.append(out_path)  #

    fp_out = open(out_path, 'w')
    s_buffer = []

    while True:
        graph = queue.get(True)
        queue_counter.decrement()
        if graph is None:
            break

        counter.decrement()

        full_index(graph, k, queue, s_buffer, haplotype,
                   prefix_sequence, queue_counter)

        for node in graph.nodes():
            s_buffer.append(gf.make_sequence(graph, node, prefix_sequence))

        graph.clear()
        graph_counter.increment()

        if len(s_buffer) > 1024:
            fp_out.write(''.join(s_buffer))
            s_buffer[:] = []

        if graph_counter.value % 1000 == 0:
            logger.info("Processed: {} subgraphs".format(graph_counter.value))

        counter.increment()

    fp_out.write(''.join(s_buffer))
    s_buffer[:] = []
    fp_out.close()

    logger.debug("Stopping: {}".format(worker_name))


def merge_files(out_file, file_list):
    with open(out_file, 'w') as outfile:
        for file in file_list:
            with open(file, 'r') as infile:
                for line in infile:
                    outfile.write(line)
            os.remove(file)


def serial_graph_wrapper(deque_collection, generator_function):
    while True:
        if len(deque_collection) == 0:
            yield next(generator_function)
        else:
            yield deque_collection.popleft()


@profile
def parallel(graph, pool, cc_queue, counter, out_path, out_file, file_list, processes, queue_counter, prefix_sequence):
    s_buffer = []

    def parallel_sleep():
        while queue_counter.value != 0 or counter.value != processes:
            sleep(0.01)

    with open(out_path, 'w') as fp_out:
        for subgraph in gf.connected_component_subgraphs(graph, s_buffer, prefix_sequence):
            cc_queue.put(subgraph)
            queue_counter.increment()
            if len(s_buffer) > 1024:
                fp_out.write(''.join(s_buffer))
                s_buffer[:] = []

        fp_out.write(''.join(s_buffer))

    parallel_sleep()

    for __ in xrange(pool._processes):
        cc_queue.put(None)
        queue_counter.increment()

    parallel_sleep()

    pool.close()
    pool.join()

    merge_files(out_file, file_list)


@profile
def serial(graph, k, out_file, haplotype, prefix_sequence=''):
    s_buffer = []
    cc_queue = deque()
    count = 0

    with open(out_file, 'w') as fp_out:
        # Iterate serially over each component in the queue
        for subgraph in serial_graph_wrapper(cc_queue, gf.connected_component_subgraphs(graph, s_buffer)):

            if len(subgraph.nodes()) == 0:
                continue

            count += 1

            if count % 1000 == 0:
                logger.info("Processed: {} subgraphs".format(count))

            full_index(subgraph, k, cc_queue, s_buffer,
                       haplotype, prefix_sequence)

            for node in subgraph.nodes():
                s_buffer.append(gf.make_sequence(subgraph, node))
            subgraph.clear()

            if len(s_buffer) > 1024:
                fp_out.write(''.join(s_buffer))
                s_buffer[:] = []

        fp_out.write(''.join(s_buffer))
        s_buffer[:] = []


def init_parallel(args):
    uniq_id = str(uuid4())
    if not os.path.exists(os.path.dirname(args.output)):
        working_directory = os.getcwd()
    else:
        working_directory = os.path.dirname(args.output)

    out_path = os.path.join(
        working_directory, "main-{}.tmp".format(uniq_id))

    cc_queue = multiprocessing.Queue()
    manager = multiprocessing.Manager()

    counter = Counter(args.processes)
    queue_counter = Counter(0)
    graph_counter = Counter(0)

    shared_prefix = manager.Value(c_char_p, args.prefix)
    file_list = manager.list()
    file_list.append(out_path)

    pool = multiprocessing.Pool(
        processes=args.processes, initializer=worker_main, initargs=((cc_queue, counter, graph_counter, queue_counter, args.kmer, file_list, working_directory, uniq_id, args.haplotype, shared_prefix)))

    return (pool, cc_queue, counter, out_path, file_list, queue_counter, shared_prefix)


@profile
def init_haplotype_graph(args):
    if args.edge_haplotype:
        graph = gf.read_gfa_edge(args.gfa)
        gf.set_to_bit_edges(graph)
    else:
        graph = gf.read_gfa(args.gfa, True)
        for node in graph.nodes():
            try:
                graph.node[node]['haplotype']
            except KeyError:
                logger.error("Node {} does not have a haplotype!".format(node))
                exit(1)
        gf.set_to_bit(graph, 2)
        gf.move_haplotype_to_edges(graph)

    return graph
