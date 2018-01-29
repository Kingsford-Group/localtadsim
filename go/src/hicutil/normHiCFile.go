package hicutil

import (
       "bufio"
       "log"
       "strconv"
       "os"
)

func ReadHiCNormFile(filename string) []float64  {

     var normvals []float64
     var x float64

     if file, err := os.Open(filename); err == nil {
     	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
	    line := scanner.Text()
	    if line == "NaN" {
	       x = 1
	    } else {
	       x, err = strconv.ParseFloat(scanner.Text(), 64)
	    }
	    normvals = append(normvals,x)
	 }

	 if err = scanner.Err(); err != nil {
	    log.Fatal(err)
	    }
	 } else {
	   log.Fatal(err)
	}
	return normvals
}