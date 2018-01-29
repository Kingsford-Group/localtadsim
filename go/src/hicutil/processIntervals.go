package hicutil

import (
//	"fmt"
//	"os"
)

func ProcessIntervals(tadlist [][]int, start int, end int) [][]int {
	// cut tadlist to just the [start, end] interval, and subtract start

	var newint [][]int
	for _,tad := range tadlist {
		if start > tad[1] { continue } // haven't gotten to start of interval yet
		if start >= tad[0]  && end >= tad[1] { // interval starts in current TAD and ends sometime afterwards
			newint = [][]int{[]int{0, tad[1]-start}}
		} else if start >= tad[0] && end < tad[1] { // interval both starts and ends in this TAD
			newint = [][]int{[]int{0,end-start}} 
		} else if start < tad[0] && end > tad[1] { // TAD fully contained in interval
			newtad := []int{tad[0] - start, tad[1] - start}
			newint = append(newint, newtad)
		} else if end >= tad[0] && end <= tad[1] { // interval ends in current TAD
			lasttad := []int{tad[0] - start, end-start}
			newint = append(newint, lasttad)
		} else if end > tad[0] { // past end of interval
			break }
		
	}
	return newint
}
