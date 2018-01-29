package hicutil

import (
	"math"
	//"fmt"
)


func CalcCondEntropy(intersects [][]int, clussizes []int, n int) float64 {

//	fmt.Println(len(clussizes))
//	fmt.Println(len(intersects))
//	fmt.Println(len(intersects[0]))
	condh := 0.0
	for _, row := range intersects {
		for j, intval := range row {
//			fmt.Println(j)
			clussize := float64(clussizes[j])
			condh += float64(intval) * math.Log(clussize)
			if intval != 0 {
				condh += - float64(intval) * math.Log(float64(intval))
			}
		}
		
	}
/*	if condh < 0 {
		fmt.Println("something wrong with cond entropy calc, negative values")
		fmt.Println(intersects)
		fmt.Println(clussizes)
		fmt.Println(n)
	}*/
	condh = condh/float64(n)
	return condh
}
