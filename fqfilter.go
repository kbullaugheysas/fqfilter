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
}

var args = Args{}

func init() {
	log.SetFlags(0)
	flag.BoolVar(&args.Invert, "invert", false, "return reads NOT in the file")
	flag.StringVar(&args.ReadsFilename, "reads", "", "filename of reads to match")
	flag.StringVar(&args.OutPrefix, "out", "", "output filename prefix (default = stdout)")

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

	// Open the inputs
	inputs := make([]AmbiReader, len(fq))
	for i, fn := range fq {
		if err := inputs[i].Open(fn); err != nil {
			log.Fatalf("Failed to open %s: %v\n", fn, err)
		}
		defer inputs[i].Close()
	}

	// Prepare the output writers
	outputs := make([]AmbiWriter, len(fq))
	for i := 0; i < len(fq); i++ {
		if args.OutPrefix == "" {
			outputs[i].Stdout()
		} else {
			fn := fmt.Sprintf("%s_%d.fq.gz", args.OutPrefix, i+1)
			if err := outputs[i].Open(fn); err != nil {
				log.Fatalf("Failed to open %s for writing: %v\n", fn, err)
			}
			defer outputs[i].Close()
		}
	}

	// Read in the list of reads
	reads := AmbiReader{}
	if err := reads.Open(args.ReadsFilename); err != nil {
		log.Fatalf("Failed to open %s: %v\n", args.ReadsFilename, err)
	}
	defer reads.Close()

	filter := make(map[string]bool)
	scanner := bufio.NewScanner(reads)
	for scanner.Scan() {
		name := scanner.Text()
		filter[name] = true
	}

	// Iterate over the inputs in sync
	inputScanners := make([]*bufio.Scanner, len(fq))
	for i := 0; i < len(fq); i++ {
		inputScanners[i] = bufio.NewScanner(inputs[i])
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
	err := func() error {
		for {
			for i := 0; i < len(fq); i++ {
				if inputScanners[i].Scan() {
					line := inputScanners[i].Text()
					if line_num%4 == 0 {
						if strings.HasPrefix(line, "@") {
							if i == 0 {
								name := line[1:len(line)]
								if args.Invert {
									_, hit := filter[name]
									enable = !hit
								} else {
									_, enable = filter[name]
								}
							}
						} else {
							return fmt.Errorf("Line %d should be a fastq header line, got: %s\n", line_num, line)
						}
					}
					if enable {
						if _, err := io.WriteString(outputs[i], line); err != nil {
							return fmt.Errorf("Failed to write line %d to output %d: %v\n", line_num, i, err)
						}
						io.WriteString(outputs[i], "\n")
						included++
					} else {
						excluded++
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
		}
	}()
	if err != nil {
		log.Fatal(err)
	}

	log.Println("included:", included)
	log.Println("excluded:", excluded)
}
