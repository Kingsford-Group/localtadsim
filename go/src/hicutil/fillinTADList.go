package hicutil

import (
	"sort"
//	"fmt"
//	"os"
)

func FillinTADList(tadlist [][]int, chrlength int) [][]int {

	var newtadlist [][]int
	if tadlist[0][0] != 0 {
		newtad := []int{0,tadlist[0][0] -1}
		newtadlist = append([][]int{newtad}, tadlist...)
	} else {
		newtadlist = tadlist
	}
	if newtadlist[len(newtadlist)-1][1] != chrlength-1 {
		newtad := []int{newtadlist[len(newtadlist)-1][1]+1, chrlength-1}
		newtadlist = append(newtadlist, newtad)
	}
	tadlistlen := len(newtadlist)
	for i:=1; i < tadlistlen; i++ {
		tadgap := newtadlist[i][0] - newtadlist[i-1][1]
		if tadgap != 1 {
			newtad := []int{newtadlist[i-1][1]+1, newtadlist[i][0]-1}
			newtadlist = append(newtadlist, newtad)
		}
	}
	sort.SliceStable(newtadlist, func(i, j int) bool {return newtadlist[i][0] < newtadlist[j][0] })
	return newtadlist

}
