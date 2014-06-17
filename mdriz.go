package main

import (
	"bufio"
	"errors"
	"flag"
	"fmt"
	"log"
	"os"
	"regexp"
	"strconv"
	"strings"
	"unsafe"
)

var linesplit = regexp.MustCompile(`\s+`)

const intsize = unsafe.Sizeof(1)

type Opts struct {
	filename *string
	depvar   *string
}

func parseArgs() *Opts {
	opts := new(Opts)
	opts.filename = flag.String("f", "", "Input file")
	opts.depvar = flag.String("d", "PHENOTYPE", "Dependency variable name")
	flag.Parse()
	return opts
}

func combinationnum(n, k int) int {
	if k > n {
		return 0
	}
	comb := make([]int, k)
	for i := range make([]struct{}, k) {
		comb[i] = i
	}
	nextComb := func(comb []int, n, k int) bool {
		i := k - 1
		comb[i]++
		for i > 0 && comb[i] == n-k+1+i {
			i--
			comb[i]++
		}
		if comb[0] > n-k {
			return false
		}
		for i = i + 1; i < k; i++ {
			comb[i] = comb[i-1] + 1
		}
		return true
	}
	count := 1
	for nextComb(comb, n, k) {
		count++
	}
	return count
}

func cleanTable(t *[2][3][3]int) {
	b := uintptr(unsafe.Pointer(t))
	for i := uintptr(0); i < 18*intsize; i += intsize {
		*(*byte)(unsafe.Pointer(b + i)) = 0
	}
}

func die(err error) {
	if err != nil {
		log.Fatal(err)
	}
}

func main() {
	opts := parseArgs()
	filer, err := os.Open(*opts.filename)
	die(err)
	defer filer.Close()
	fscanner := bufio.NewScanner(filer)
	samplesize := 0
	for fscanner.Scan() {
		samplesize++
	}
	samplesize -= 1
	filer.Seek(0, 0)
	fscanner = bufio.NewScanner(filer)
	fscanner.Scan()
	varnames := linesplit.Split(fscanner.Text(), -1)
	varnum := len(varnames)
	if varnum < 2 {
		log.Fatal(errors.New("Can not read variable names on first line"))
	}

	var depname string
	var depindex int
	for i, dep := range varnames {
		if dep == *opts.depvar {
			depname = *opts.depvar
			depindex = i
		}
	}
	if depname == "" {
		log.Fatal(errors.New("Can not find dependency variable name on first line"))
	}

	datamat := make([][]byte, 0, varnum)
	for _ = range make([]struct{}, varnum) {
		datamat = append(datamat, make([]byte, samplesize))
	}
	{
		c := 0
		for fscanner.Scan() {
			gtStrings := linesplit.Split(fscanner.Text(), -1)
			for i, j := range gtStrings {
				gt, err := strconv.ParseInt(j, 10, 8)
				die(err)
				datamat[i][c] = byte(gt)
			}
			c++
		}
	}
	var casePortion float64
	{
		var sum int
		for _, j := range datamat[depindex] {
			sum += int(j)
		}
		casePortion = float64(sum) / float64(samplesize-sum)
	}

	combnum := combinationnum(varnum, 2)
	//classTable := make(map[[2]int]*[3][3]bool, combnum)
	output := make([][]byte, 0, combnum)
	outputHeader := make([]string, 0, combnum)

	output = append(output, datamat[depindex])
	outputHeader = append(outputHeader, depname)
	{
		classTable := [3][3]bool{}
		table := [2][3][3]int{}
		for i1 := 0; i1 < varnum-1; i1++ {
			if i1 == depindex || i1+1 == depindex {
				continue
			}
			for i2 := i1 + 1; i2 < varnum; i2++ {
				//	key := [2]int{i1, i2}
				cleanTable(&table)
				for j := range make([]struct{}, samplesize) {
					table[datamat[depindex][j]][datamat[i1][j]][datamat[i2][j]]++
				}
				//classTable[key] = &[3][3]bool{}
				for a1 := 0; a1 < 3; a1++ {
					for a2 := 0; a2 < 3; a2++ {
						if float64(table[1][a1][a2]) >= casePortion*float64(table[0][a1][a2]) {
							//classTable[key][a1][a2] = true
							classTable[a1][a2] = true

						} else {
							classTable[a1][a2] = false
						}
					}
				}
				newvar := make([]byte, samplesize)
				for k := 0; k < samplesize; k++ {
					if classTable[datamat[i1][k]][datamat[i2][k]] {
						newvar[k] = 1
					} else {
						newvar[k] = 0
					}
				}
				output = append(output, newvar)
				outputHeader = append(outputHeader, fmt.Sprintf("%s_x_%s", varnames[i1], varnames[i2]))
			}
		}
	}
	os.Stdout.WriteString(strings.Join(outputHeader, "\t"))
	os.Stdout.WriteString("\n")

	line := make([]string, len(output))
	for k := 0; k < samplesize; k++ {
		for m := 0; m < len(output); m++ {
			line[m] = fmt.Sprintf("%d", output[m][k])
		}
		os.Stdout.WriteString(strings.Join(line, "\t"))
		os.Stdout.WriteString("\n")
	}
}
