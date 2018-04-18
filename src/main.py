from chop import *


def main():
    logging.basicConfig(stream=sys.stdout,
                        format='%(asctime)s - %(levelname)s:%(name)s:%(message)s')

    parser = argparse.ArgumentParser(
        description='Index k-length paths in provided graph in GFA format with or without haplotype constraints.')
    parser.add_argument('-g', '--gfa', dest='gfa',
                        help='Input GFA file', type=str, required=True)
    parser.add_argument('-k', dest='kmer',
                        help='k-mer value', type=int, required=True)
    parser.add_argument('-o', dest='output',
                        help='Output file', type=str, required=True)
    parser.add_argument('-p', dest='processes',
                        help='Number of parallel processes to start', type=int)
    parser.add_argument('-s', '--prefix', dest='prefix',
                        help='Optional prefix for each path entry', type=str, default='')
    parser.add_argument('-H', '--haplotype', dest='haplotype',
                        help='Use haplotype information', action='store_true')
    parser.add_argument('-e', '--edge', dest='edge_haplotype',
                        help='Haplotypes are already encoded on the edges', action='store_true')
    parser.add_argument('-l', dest='log_level', help='Set the logging level',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'])
    parser.add_argument('--version', action='version', version='%(prog)s: 0.2')
    args = parser.parse_args()

    if args.log_level:
        logging.getLogger().setLevel(getattr(logging, args.log_level))

    if args.processes and args.processes > 1:
        (pool, cc_queue, counter, out_path, file_list,
         queue_counter, shared_prefix) = init_parallel(args)

    if args.haplotype:
        graph = init_haplotype_graph(args)
    else:
        graph = gf.read_gfa(args.gfa, False)

    gf.set_interval_nodes(graph)
    graph.k_extend = args.kmer - 1
    graph.max_id = max([i for i in graph.nodes()])

    first_sweep(graph)

    if args.processes and args.processes > 1:
        parallel(graph, pool, cc_queue, counter, out_path, args.output,
                 file_list, args.processes, queue_counter, shared_prefix.value)
    else:
        serial(graph, args.kmer, args.output, args.haplotype, args.prefix)


if __name__ == '__main__':
    main()
