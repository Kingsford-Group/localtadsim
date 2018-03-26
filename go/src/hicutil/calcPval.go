package hicutil

import (
	"math"
	"math/rand"
//	"sync"
//	"runtime"
//	"fmt"
//	"os"
)

func CalcPval(intvl1 [][]int, intvl2 [][]int, n int, vival float64, nshuffles int) float64 {
// more efficient to calculate VI as entropy(intvl1) + entropy(intvl2) - 2*mutinfo, because only need to recalculate mutual info on each iteration
	shuffvi1 := make([]float64, nshuffles)
	shuffvi2 := make([]float64, nshuffles)

	h1 := CalcEntropy(intvl1,n)
	h2 := CalcEntropy(intvl2,n)
	//n := intvl1[len(intvl1)-1][1] - intvl1[0][0]+1

	clus1sizes := make([]int, len(intvl1))
        for c,clus := range intvl1 {
		clus1sizes[c] = clus[1]-clus[0]+1 }
	clus2sizes := make([]int, len(intvl2))
	for c,clus := range intvl2 {
		clus2sizes[c] = clus[1]-clus[0]+1 }

	// parallelize/make these concurrent??
//	var wg sync.WaitGroup
//	wg.Add(nshuffles)
//	runtime.GOMAXPROCS(nshuffles)
	for i := 0; i < nshuffles; i++ {
//		go func (i int) {
//		defer wg.Done()

		// randomly shuffle domain lengths in each list
		newlist1,permsizes1 := shuffledoms(clus1sizes)
		newlist2,permsizes2 := shuffledoms(clus2sizes)
		
		//calc VI for newly shuffled domains
		overlaps1 := CalcOverlaps(intvl1, newlist2)
		overlaps2 := CalcOverlaps(newlist1, intvl2)
		mutinfo1 := CalcMutInfo(overlaps1, clus1sizes, permsizes2, n)
		mutinfo2 := CalcMutInfo(overlaps2, permsizes1, clus2sizes, n)
		
		shuffvi1[i] = (h1+h2-2*mutinfo1)/math.Log(float64(n)) // divide by log(n) to normalize
		shuffvi2[i] = (h1+h2-2*mutinfo2)/math.Log(float64(n)) // divide by log(n) to normalize
		
			/*// test VI calculation
		condhvi1 := CalcCondEntropy(transpose(overlaps1), clus1sizes, n) + CalcCondEntropy(overlaps1, permsizes2, n)
		condhvi2 := CalcCondEntropy(transpose(overlaps2), permsizes1, n) + CalcCondEntropy(overlaps2, clus2sizes, n)
		
		if math.Abs(condhvi1 - shuffvi1[i]) > 1e-10 || math.Abs(condhvi2 - shuffvi2[i]) > 1e-10 {
			fmt.Println("VI calculation is wrong")
			fmt.Println(shuffvi1[i], condhvi1)
			fmt.Println(shuffvi2[i], condhvi2)
			fmt.Println(intvl1, newlist2)
			fmt.Println(intvl2, newlist1)
			os.Exit(1)
		}*/

		//fmt.Println("origvi =",vival,"newvi1 =",shuffvi1[i],"newvi2 =",shuffvi2[i])
//		}(i)
	}
//	wg.Wait()
	// find how many times shuffvi values are less than given vival
	count1 := 0
	count2 := 0
	for i:=0; i< nshuffles; i++ {
		if shuffvi1[i] - vival < 1e-10 { count1++ }
		if shuffvi2[i] - vival < 1e-10 { count2++ }
	}
	pval := (float64(count1+1)/float64(nshuffles+1) + float64(count2+1)/float64(nshuffles+1))/2.0
/*	if pval < 0.05 && (len(intvl1) == 1 || len(intvl2) == 1) {
		fmt.Println(intvl1)
		fmt.Println(intvl2)
		fmt.Println(pval, count1, count2)
		fmt.Println(n,vival,nshuffles)
		fmt.Println(shuffvi1)
		fmt.Println(shuffvi2)
		//fmt.Println(querypt.start, querypt.end)
                os.Exit(1)
        }*/
	return pval
}

func shuffledoms(clussizes []int) ([][]int, []int) {

	permsizes := make([]int, len(clussizes))
	perm := rand.Perm(len(clussizes))
	for j,v := range perm {
		permsizes[v] = clussizes[j]
	}
	// turn shuffled lists of lengths back into TAD lists
	newlist := make([][]int, len(clussizes))
	newlist[0] = []int{0,permsizes[0]-1}
	for j := 1; j < len(permsizes); j++ {
		tadstart := newlist[j-1][1]+1
		newlist[j] = []int{tadstart, tadstart+permsizes[j]-1}
	}

	return newlist, permsizes
}
