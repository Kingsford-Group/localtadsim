package hicutil

import (
	"path/filepath"
	"log"
	"math"
	"strings"
	"strconv"
	"sort"
//	"fmt"
)

func ChooseGamma(tadlen float64,fileseed string, res int) ([][]int, float64) {

	filelist,err := filepath.Glob(fileseed+"*")
//	fmt.Println(filelist)
	if err != nil {
		log.Fatal(err)
	}
	tadsets := make(map[float64][][]int)
	for _,file := range filelist {
		fparts := strings.Split(file, ".")
		gammastart := findIndex(fparts, "gamma") + 1
		checkconsensus := findIndex(fparts, "consensus")
		if gammastart == 0 && checkconsensus == -1 {
			log.Println("Armatus file doesn't contain the word 'gamma', make sure file path is correct")
			log.Fatal(err)
		} else if checkconsensus > -1 {
			continue
		}
		gammaslice := fparts[gammastart:len(fparts)-2]
		gammastr := strings.Join(gammaslice,".")
		gamma,err := strconv.ParseFloat(gammastr, 64)
		if err != nil {
			log.Fatal(err)
		}
		tadsets[gamma] = ReadTADFile(file, res)
	}
	var disttotadlen,optgamma float64
	for g,tadlist := range tadsets {
//		var sumtadsize int
		var tadsizes []int
		for _,tad := range tadlist {
			//put TAD lengths into slice
			tadsizes = append(tadsizes, tad[1]-tad[0])
//			sumtadsize += tad[1]-tad[0]
		}
		// find median of slice
		medtadlen := median(tadsizes)
//		meantadlen := float64(sumtadsize)/float64(len(tadlist))
//		currdist := math.Abs(meantadlen - tadlen)
		currdist := math.Abs(medtadlen - tadlen)
		if (disttotadlen == 0 && medtadlen != tadlen) || currdist < disttotadlen {
			disttotadlen = currdist
			optgamma = g
		}
	}
//	fmt.Println(fileseed,optgamma)
	return tadsets[optgamma], optgamma
}


func findIndex(s []string, exp string) int {
	for idx,v := range s {
		if v == exp {
			return idx
		}
	}
	return -1
}

func median(a []int) float64 {

	sort.Ints(a)
	midpt := len(a)/2
	if math.Mod(float64(len(a)),2.0) == 0 {
		return float64(a[midpt] + a[midpt-1])/2.0
	} else {
		return float64(a[midpt])
	}
}
