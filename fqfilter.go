package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"strings"
)

/* This program takes on one or two (in the case of paried end data) fq files
 * and returns a subset of the reads */

type Args struct {
	Invert        bool
	ReadsFilename string
	OutPrefix     string
	Limit         int
	Tab           bool
	ShortName     bool
}

var args = Args{}

func init() {
	log.SetFlags(0)
	flag.BoolVar(&args.Invert, "invert", false, "return reads NOT in the file")
	flag.BoolVar(&args.Tab, "tab", false, "print sequence as tabular output (readName, read1, read2)")
	flag.BoolVar(&args.ShortName, "short-name", false, "use just the first space-separated word of the read name")
	flag.StringVar(&args.ReadsFilename, "reads", "", "filename of reads to match")
	flag.StringVar(&args.OutPrefix, "out", "", "output filename prefix (default = stdout)")
	flag.IntVar(&args.Limit, "limit", 0, "output only the first LIMIT matches")

	flag.Usage = func() {
		log.Println("usage: fqfilter [options] unaligned_1.fq.gz unaligned_2.fq.gz")
		flag.PrintDefaults()
	}
}

/* Provide an ambidexterous interface to files to read that may be gzipped */
type AmbiReader struct {
	fp *os.File
	gz *gzip.Reader
	r  io.Reader
}

func (a AmbiReader) Read(b []byte) (n int, err error) {
	return a.r.Read(b)
}

func (a *AmbiReader) Open(fn string) error {
	if a.r != nil {
		return fmt.Errorf("AmbiReader already open")
	}
	var err error
	// If no filename is given, then read from stdin
	if fn == "" {
		a.r = os.Stdin
		return nil
	}
	a.fp, err = os.Open(fn)
	if err != nil {
		return err
	}
	if strings.HasSuffix(fn, ".gz") {
		a.gz, err = gzip.NewReader(a.fp)
		if err != nil {
			return err
		}
		a.r = a.gz
	} else {
		a.r = a.fp
	}
	return nil
}

func (a *AmbiReader) Close() error {
	if a.gz != nil {
		if err := a.gz.Close(); err != nil {
			return err
		}
	}
	if err := a.fp.Close(); err != nil {
		return err
	}
	return nil
}

/* Provide an ambidexterous interface to files to write that may be gzipped */
type AmbiWriter struct {
	fp *os.File
	gz *gzip.Writer
	r  io.Writer
}

func (a AmbiWriter) Write(b []byte) (n int, err error) {
	return a.r.Write(b)
}

func (a *AmbiWriter) Close() error {
	if a.gz != nil {
		if err := a.gz.Close(); err != nil {
			return err
		}
	}
	if err := a.fp.Close(); err != nil {
		return err
	}
	return nil
}

func (a *AmbiWriter) Open(fn string) error {
	if a.r != nil {
		return fmt.Errorf("AmbiWriter already open")
	}
	var err error
	// If no filename is given, then read from stdin
	if fn == "" {
		a.r = os.Stdout
		return nil
	}
	a.fp, err = os.Create(fn)
	if err != nil {
		return err
	}
	if strings.HasSuffix(fn, ".gz") {
		a.gz = gzip.NewWriter(a.fp)
		a.r = a.gz
	} else {
		a.r = a.fp
	}
	return nil
}

func (a *AmbiWriter) Stdout() {
	a.r = os.Stdout
}

func main() {
	flag.Parse()
	fq := flag.Args()

	if args.ReadsFilename == "" {
		log.Fatal("Must provide -reads <file> argument")
	}

	if len(fq) == 0 {
		log.Fatal("Must specify at least one fastq file")
	}

	// Open the inputs
	inputs := make([]AmbiReader, len(fq))
	for i, fn := range fq {
		if err := inputs[i].Open(fn); err != nil {
			log.Fatalf("Failed to open %s: %v\n", fn, err)
		}
		defer inputs[i].Close()
	}

	var outputs []AmbiWriter

	if args.Tab {
		if args.OutPrefix != "" {
			log.Fatal("Tabular output only supports writing to stdout")
		}
	} else {
		// Prepare the output writers

		outputs = make([]AmbiWriter, len(fq))
		for i := 0; i < len(fq); i++ {
			if args.OutPrefix == "" {
				outputs[i].Stdout()
			} else {
				var fn string
				if len(fq) == 1 {
					fn = fmt.Sprintf("%s.fq.gz", args.OutPrefix)
				} else {
					fn = fmt.Sprintf("%s_%d.fq.gz", args.OutPrefix, i+1)
				}
				if err := outputs[i].Open(fn); err != nil {
					log.Fatalf("Failed to open %s for writing: %v\n", fn, err)
				}
				defer outputs[i].Close()
			}
		}
	}

	// Read in the list of reads
	reads := AmbiReader{}
	readsFn := args.ReadsFilename
	if readsFn == "stdin" {
		readsFn = ""
	}
	if err := reads.Open(readsFn); err != nil {
		log.Fatalf("Failed to open %s: %v\n", args.ReadsFilename, err)
	}
	defer reads.Close()

	filter := make(map[string]bool)
	scanner := bufio.NewScanner(reads)
	for scanner.Scan() {
		name := scanner.Text()
		if args.ShortName {
			name = strings.Fields(name)[0]
		}
		filter[name] = true
	}

	// Iterate over the inputs in sync
	inputScanners := make([]*bufio.Scanner, len(fq))
	for i := 0; i < len(fq); i++ {
		inputScanners[i] = bufio.NewScanner(inputs[i])
		/* Make sure we have a large buffer for long sequences */
		buf := make([]byte, 0, 1024*1024)
		inputScanners[i].Buffer(buf, 10*1024*1024)
	}
	var enable bool
	if args.Invert {
		enable = true
	} else {
		enable = false
	}
	line_num := 0
	included := 0
	excluded := 0
	var name string
	sequences := make([]string, len(fq))
	err := func() error {
		for {
			for i := 0; i < len(fq); i++ {
				if inputScanners[i].Scan() {
					line := inputScanners[i].Text()
					if line_num%4 == 0 {
						if strings.HasPrefix(line, "@") {
							if i == 0 {
								name = line[1:len(line)]
								if args.ShortName {
									name = strings.Fields(name)[0]
								}
								_, enable = filter[name]
								if args.Invert {
									enable = !enable
								}
							}
						} else {
							return fmt.Errorf("Line %d should be a header line, got: %s\n", line_num, line)
						}
					}
					if line_num%4 == 1 {
						sequences[i] = line
						if i == 0 {
							if enable {
								included++
							} else {
								excluded++
							}
						}
					}
					if enable {
						if args.Tab {
							if i+1 == len(fq) && line_num%4 == 1 {
								outputLine := name
								for j := 0; j < len(fq); j++ {
									outputLine = outputLine + "\t" + sequences[j]
								}
								fmt.Println(outputLine)
							}
						} else {
							if _, err := io.WriteString(outputs[i], line+"\n"); err != nil {
								return fmt.Errorf("Failed to write line %d to output %d: %v\n", line_num, i, err)
							}
						}
					}
				} else {
					if i == 0 {
						return nil
					} else {
						return fmt.Errorf("Expecting scanner %d to be able to scan\n", i)
					}
				}
			}
			line_num++
			if args.Limit > 0 && included >= args.Limit {
				log.Println("reached limit")
				return nil
			}
		}
	}()
	if err != nil {
		log.Fatal(err)
	}

	log.Println("included:", included)
	log.Println("excluded:", excluded)
}
