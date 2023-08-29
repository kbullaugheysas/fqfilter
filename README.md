# fqfilter
Utility for filtering a fastq file based on a list of read names.

    usage: fqfilter [options] unaligned_1.fq.gz unaligned_2.fq.gz
      -invert
            return reads NOT in the file
      -limit int
            output only the first LIMIT matches
      -out string
            output filename prefix (default = stdout)
      -reads string
            filename of reads to match
      -short-name
            use just the first space-separated word of the read name
      -tab
            print sequence as tabular output (readName, read1, read2)
