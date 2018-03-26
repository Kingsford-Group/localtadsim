
package hicutil

import(
	"bufio"
	"log"
	"strings"
	"strconv"
	"os"
	"math"
//	"fmt"
)


func ReadTADFile(filename string, res int) [][]int  {

	var tadlist [][]int
	var start,end int

	if file, err := os.Open(filename); err == nil {
		defer file.Close()

		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			tad := make([]int, 2)
			line := strings.Fields(scanner.Text())
			start,err = strconv.Atoi(line[1])
			start = start/res
			end,err = strconv.Atoi(line[2])
			// (end+1)/res -1 for Armatus - for others just do end/res-1
			if math.Mod(float64(end), float64(res)) != 0 {
				end = (end+1)/res - 1
			} else {
				end = end/res - 1
			}	
			tad[0] = start
			tad[1] = end
			tadlist = append(tadlist,tad)

		}
		
		if err = scanner.Err(); err != nil {
			log.Fatal(err)
		}
		} else {
			log.Fatal(err)
		}
	if tadlist[0][0] > tadlist[len(tadlist)-1][0] {
                // TADs are listed backwards, need to flip 
                tadlist = flipTADs(tadlist)
        }
	
	return tadlist
}

func flipTADs(tadlist [][]int) [][]int {

        flippedtads := make([][]int, len(tadlist))
        for i,tad := range tadlist {
                flippedtads[len(tadlist)-1-i] = tad
        }
        return flippedtads
}

