package hicutil

import (
	"fmt"
)

func CalcOverlaps(cluslist1 [][]int, cluslist2 [][]int) [][]int {
// inputs should just be 2 lists of intervals, both starting at 0 and ending at the same n
// output is matrix of overlaps, where [i][j] entry is the size of the intersection of cluster i from cluslist1 with cluster j from cluslist2

	overlaps := make([][]int, len(cluslist1))
	for i,clus1 := range cluslist1 {
		overlaps[i] = make([]int, len(cluslist2))
		for j,clus2 := range(cluslist2) {
			if clus2[0] > clus1[1] || clus1[0] > clus2[1] {
				overlaps[i][j] = 0 // no overlap
			} else if clus1[0] <= clus2[0] && clus2[0] <= clus1[1] && clus2[1] > clus1[1] {
				overlaps[i][j] = clus1[1] - clus2[0] + 1 //partial overlap (clus1 starts first)
			} else if clus2[0] <= clus1[0] && clus1[0] <= clus2[1] && clus1[1] > clus2[1] {
				overlaps[i][j] = clus2[1] - clus1[0] + 1 //partial overlap (clus2 starts first)
			} else if clus2[0] >= clus1[0] && clus2[1] <= clus1[1] {
				overlaps[i][j] = clus2[1] - clus2[0] + 1 //clus2 fully contained in clus1
			} else if clus1[0] >= clus2[0] && clus1[1] <= clus2[1] {
				overlaps[i][j] = clus1[1] - clus1[0] + 1 //clus1 fully contained in clus2
			} else {
				fmt.Println("this wasn't supposed to happen")
			}
		}
	}
	return overlaps
}


func CalcOverlapsPtr(cluslist1 [][]int, cluslist2 [][]int, overlaps *[][]int) {
// inputs should just be 2 lists of intervals, both starting at 0 and ending at the same n
// output is matrix of (*overlaps), where [i][j] entry is the size of the intersection of cluster i from cluslist1 with cluster j from cluslist2

	for i,clus1 := range cluslist1 {
		for j,clus2 := range(cluslist2) {
			if clus2[0] > clus1[1] || clus1[0] > clus2[1] {
				(*overlaps)[i][j] = 0 // no overlap
			} else if clus1[0] <= clus2[0] && clus2[0] <= clus1[1] && clus2[1] > clus1[1] {
				(*overlaps)[i][j] = clus1[1] - clus2[0] + 1 //partial overlap (clus1 starts first)
			} else if clus2[0] <= clus1[0] && clus1[0] <= clus2[1] && clus1[1] > clus2[1] {
				(*overlaps)[i][j] = clus2[1] - clus1[0] + 1 //partial overlap (clus2 starts first)
			} else if clus2[0] >= clus1[0] && clus2[1] <= clus1[1] {
				(*overlaps)[i][j] = clus2[1] - clus2[0] + 1 //clus2 fully contained in clus1
			} else if clus1[0] >= clus2[0] && clus1[1] <= clus2[1] {
				(*overlaps)[i][j] = clus1[1] - clus1[0] + 1 //clus1 fully contained in clus2
			} else {
				fmt.Println("this wasn't supposed to happen")
			}
		}
	}
}
