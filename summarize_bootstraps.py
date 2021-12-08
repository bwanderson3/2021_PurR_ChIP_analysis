if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='commands', dest='command')

    # parse
    parse_parser = subparsers.add_parser('parse', help="create a sampling\
            object from a sam file")
    parse_parser.add_argument('samfile', help="Input samfile, this tool does no\
            filtering and will consider every line in the file. Accepts input\
            from stdin if '-' is specified here.")
    parse_parser.add_argument('outpre', help="output prefix to np.save the\
            sampling object data that is created")
    parse_parser.add_argument('--paired', action="store_true", help="Consider\
            the sam file as paired. If this flag is specified then the sam file\
            MUST be pre-filtered to have only ONE alignment per pair. Further,\
            there must be NO unpaired reads in the file and the reads must be\
            sorted by read name.")

    # sample
    sample_parser = subparsers.add_parser('sample', help="sample coverage from a\
            sampling object.")
    sample_parser.add_argument('samplerfile', help="output file from using parse\
            on the samfile of interest")
    sample_parser.add_argument('outpre', help="output file to np.save the\
            numpy array that is created")
    sample_parser.add_argument('--reference_genome', type=str, help="fasta\
        file containing reference genome. Will be used to infer genome size")
    # sample_parser.add_argument('array_size',type=int, help="length of genome")
    sample_parser.add_argument('--num_samples', type=int, default=1,
    help="number of full samples to pull from the sampler, default is 1")
    sample_parser.add_argument('--num_reads', type=int, default=None,
    help="number of reads to pull for each sample. Default is the size of\
            sampling object.")
    sample_parser.add_argument('--identity', action="store_true",
            help="write an array of the actual sample without sampling, ignores\
                    num_reads, num_samples and seed options")
    sample_parser.add_argument('--resolution', type=int, default=1,
            help="only report every x basepairs, default=1")
    sample_parser.add_argument('--seed', type=int, default=1234,
            help="psuedo-random number generator seed, default=1234")
    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s %(message)s',level=logging.INFO)
    if args.command == "parse":
        logging.info("Preparing sampler object {}".format(args.outpre + ".npy"))
        if args.samfile == "-":
            logging.info("Reading stream from stdin")
            f = sys.stdin
        else:
            logging.info("Reading file {}".format(args.samfile))
            f = open(args.samfile, mode="r")
        if args.paired:
            sampler = create_read_list_paired(f)
        else:
            sampler = create_read_list(f)
        f.close()
        sampler.sort_reads()
        sampler.save_data(args.outpre)

    elif args.command == "sample":
        array_size = len(SeqIO.read(args.reference_genome, 'fasta'))

        logging.info("Sampling from file {}".format(args.outpre + ".npy"))
        logging.info("Using seed {}".format(args.seed))
        prng = np.random.RandomState(args.seed)
        array = np.zeros((int(floor(array_size/args.resolution)), args.num_samples), dtype=np.int32)
        sampler = ReadSampler()
        sampler.load_data(args.samplerfile)
        if args.identity:
            for read in sampler.reads:
                map_read(array, read, args.resolution)
            np.save(args.outpre, array)
        else:
            for i in range(args.num_samples):
                begin = time.time()
                if args.num_reads:
                    num_reads = args.num_reads
                else:
                    num_reads = sampler.total
                # array = sample(sampler, num_reads, array[:,i], args.resolution, prng)
                sample(sampler, num_reads, array[:,i], args.resolution, prng)
                finish = time.time()
                logging.info("Sample {} took {} seconds".format(i, finish-begin))
            np.save(args.outpre, array)
