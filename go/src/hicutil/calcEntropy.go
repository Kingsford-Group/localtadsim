package hicutil

import (
	"math"
	"fmt"
//	"os"
)

func CalcEntropy(clusters [][]int, intn int) float64 {

	n := float64(intn)
	//n := float64(clusters[len(clusters)-1][1] - clusters[0][0] + 1)
	entropy := 0.0
	for _,clus := range clusters {
		clussize := float64(clus[1]-clus[0]+1)
		entropy += -(clussize/n)*math.Log(clussize/n)
	}
	if math.IsNaN(entropy) {
		fmt.Println("entropy is NaN")
		fmt.Println(clusters)
	}
	return entropy
}
