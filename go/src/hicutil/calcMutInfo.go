package hicutil

import (
	"math"
	"fmt"
//	"os"
)

func CalcMutInfo(overlaps [][]int, sizes1 []int, sizes2 []int, intn int) float64 {

	mutinfo := 0.0
	n := float64(intn)
	for i,row := range overlaps {
		for j,val := range row {
			jointpr := float64(val)/n
			px := float64(sizes1[i])/n
			py := float64(sizes2[j])/n
			if jointpr != 0 {
				mutinfo += jointpr*math.Log(jointpr/(px*py))
			}
		}
	}
	if math.IsNaN(mutinfo) || mutinfo < 0 {
		fmt.Println("mutinfo is NaN or negative")
		fmt.Println("mutinfo =",mutinfo)
		fmt.Println(overlaps)
		fmt.Println(sizes1)
		fmt.Println(sizes2)
	}
	return mutinfo
}
