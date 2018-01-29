package main

import (
	"math"
	"hicutil"
//	"fmt"
)

func CorrectTADoverlap( tadlists [][][]int, i int, start int, endtad []int, dpmap map[int]map[int]condentropy ) float64 {
	// there is a TAD that overlaps the boundary of endtad, so have to correct conditional entropy term
	n := endtad[1] - start + 1
	var noti int
	if i == 0 {noti = 1} else {noti = 0}
	// find the overlapping TAD (first TAD containing the beginning of endtad)
	var ovrlpi int
	for idx,tads := range tadlists[noti] {
		if endtad[0] <= tads[1] && endtad[0] >= tads[0] {
			ovrlpi = idx
			break
		}
	}
	overlaptad := tadlists[noti][ovrlpi]
	/*fmt.Println("endtad =",endtad)
	fmt.Println("overlaptad =",overlaptad)
	fmt.Println("prev tad =",tadlists[noti][ovrlpi-1])*/
	otad1 := []int{max(overlaptad[0], start), endtad[0]-1}
	//fmt.Println(otad1)
	n1 := otad1[1] - otad1[0]+1
	otad2 := []int{endtad[0], min(overlaptad[1], endtad[1])}
//	fmt.Println(otad2)
	n2 := otad2[1] - otad2[0]+1
//	otadn := otad2[1]-otad1[0]+1
	scaledint1 := hicutil.ProcessIntervals(tadlists[i], otad1[0], otad1[1])
	scaledint2 := hicutil.ProcessIntervals(tadlists[i], otad2[0], otad2[1])
	// calculate intersects and terms to remove
	intersects1 := hicutil.CalcOverlaps([][]int{[]int{0,otad1[1]-otad1[0]}}, scaledint1)
	intersects2 := hicutil.CalcOverlaps([][]int{[]int{0,otad2[1]-otad2[0]}}, scaledint2)
	sec1remove := 0.0
	for _,ival := range intersects1[0] {
		if ival != 0 {
			sec1remove += float64(ival)*math.Log(float64(n1)) - float64(ival)*math.Log(float64(ival))
		}
	}
	sec2remove := 0.0
	for _,ival := range intersects2 {
		if ival[0] != 0 {
			sec2remove += float64(ival[0])*math.Log(float64(n2)) - float64(ival[0])*math.Log(float64(ival[0]))
		}
	}
	var newcondh float64
	if i == 0 {
		newcondh = (float64(endtad[0] - start)/float64(n))*dpmap[start][endtad[0]-1].condh2 + (float64(endtad[1] - endtad[0]+1)/float64(n))*dpmap[endtad[0]][endtad[1]].condh2 
	} else {
		newcondh = (float64(endtad[0] - start)/float64(n))*dpmap[start][endtad[0]-1].condh1 + (float64(endtad[1] - endtad[0]+1)/float64(n))*dpmap[endtad[0]][endtad[1]].condh1
	}
	otadn := otad2[1] - otad1[0] + 1
	otadint := hicutil.ProcessIntervals(tadlists[i], otad1[0], otad2[1])
	newintersects := hicutil.CalcOverlaps(otadint, [][]int{[]int{0, otadn}})
	newterm := hicutil.CalcCondEntropy(newintersects, []int{otadn}, otadn)
	newcondh2 := newcondh -(1/float64(n))*sec1remove - (1/float64(n))*sec2remove + (float64(otadn)/float64(n))*newterm
	return newcondh2
}



func max(a int, b int) int {
	if a >= b {
		return a
	} else {
		return b
	}
}

func min(a int, b int) int {
	if a <= b {
		return a 
	} else {
		return b
	}
}
