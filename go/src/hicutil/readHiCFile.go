package hicutil

import (
	"bufio"
	"log"
	"strings"
	"strconv"
	"os"
	"math"
)

// define a pair (tuple) type for keys to Hi-C maps
type Pair struct {
     A,B interface{}
}

func ReadHiCFile(filename string, res int, normvals []float64) (map[Pair]float64,int) {

     hicmap := make(map[Pair]float64)
     var row,col int
     var val float64
     var norm bool
     var chrlength int

     if len(normvals) > 0 {
     	norm = true
	} else {
       norm = false
       }

     if file, err := os.Open(filename); err == nil {
     	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
	    vals := strings.Fields(scanner.Text())
	    row, err = strconv.Atoi(vals[0])
	    row = row/res
	    col, err = strconv.Atoi(vals[1])
	    col = col/res
	    
	    if row > chrlength {
	       chrlength = row
	    }
	    if col > chrlength {
	       chrlength = col
	    }

	    val, err = strconv.ParseFloat(vals[2],64)
	    if norm {
	       val = float64(val) / (float64(normvals[row])*float64(normvals[col]))
	    }
	hicmap[Pair{row,col}] = math.Log(val)
//	hicmap[Pair{row,col}] = val
	}

	if err = scanner.Err(); err != nil {
	   log.Fatal(err)
	   }
	} else {
	  log.Fatal(err)
	}
	return hicmap, chrlength
}
