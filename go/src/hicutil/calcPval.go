package hicutil

import (
	"math"
	"math/rand"
//	"sync"
//	"runtime"
//	"fmt"
//	"os"
)

func CalcPval(intvl1 [][]int, intvl2 [][]int, n int, vival float64, convcond float64, r *rand.Rand) float64 {
// more efficient to calculate VI as entropy(intvl1) + entropy(intvl2) - 2*mutinfo, because only need to recalculate mutual info on each iteration
	//shuffvi1 := make([]float64, nshuffles)
	//shuffvi2 := make([]float64, nshuffles)

	//fmt.Println("convcond =",convcond)

	//var shuffvi []float64
	var pval1 float64
	var pval2 float64	

	h1 := CalcEntropy(intvl1,n)
	h2 := CalcEntropy(intvl2,n)
	//n := intvl1[len(intvl1)-1][1] - intvl1[0][0]+1

	clus1sizes := make([]int, len(intvl1))
        for c,clus := range intvl1 {
		clus1sizes[c] = clus[1]-clus[0]+1 }
	clus2sizes := make([]int, len(intvl2))
	for c,clus := range intvl2 {
		clus2sizes[c] = clus[1]-clus[0]+1 }
	
	overlaps1 := make([][]int, len(intvl1))
	for i,_ := range intvl1 {
		overlaps1[i] = make([]int, len(intvl2))
		}
	overlaps2 := make([][]int, len(intvl1))
	for i,_ := range intvl1 {
                overlaps2[i] = make([]int, len(intvl2))
        }

	count := 0
	shuffnum := 0
	prevpvaldiffs := []float64{1000.0,1000.0,1000.0,1000.0}
	prevpval := 1000.0
	keepshuffling := true
 
	for keepshuffling {
		shuffnum++

		// randomly shuffle domain lengths in each list
		newlist,permsizes := shuffledoms(clus1sizes, r)
		
		//calc VI for newly shuffled domains
		CalcOverlapsPtr(newlist, intvl2, &overlaps1) //change
		mutinfo := CalcMutInfo(overlaps1, permsizes, clus2sizes, n)
		
		shuffvi := (h1+h2-2*mutinfo)/math.Log(float64(n)) // divide by log(n) to normalize
		
		// calc current pval, and pvaldiffs
		if shuffvi - vival < 1e-10 { count++ }
		pval1 = float64(count+1)/float64(shuffnum+1)

		prevpvaldiffs = append(prevpvaldiffs[1:],[]float64{math.Abs(prevpval-pval1)}...)
		// if last 5 p-values are within convergence condition, we're done shuffling
		keepshuffling = false
		for _,pvaldiff := range prevpvaldiffs {
			if pvaldiff > convcond { keepshuffling = true }
		}
		prevpval = pval1
	}

	// re-initialize variables
	shuffnum = 0
	count = 0
	prevpvaldiffs = []float64{1000.0,1000.0,1000.0,1000.0}
	prevpval = 1000.0
	keepshuffling = true

	for keepshuffling {
		shuffnum++

		newlist,permsizes := shuffledoms(clus2sizes, r)
		CalcOverlapsPtr(intvl1, newlist, &overlaps2)
		mutinfo := CalcMutInfo(overlaps2, clus1sizes, permsizes, n)
		
		shuffvi := (h1+h2-2*mutinfo)/math.Log(float64(n)) // divide by log(n) to normalize
		
		// calc current pval, and pvaldiff
		if shuffvi - vival < 1e-10 { count++ }
		pval2 = float64(count+1)/float64(shuffnum+1)

		prevpvaldiffs = append(prevpvaldiffs[1:],[]float64{math.Abs(prevpval-pval2)}...)
		keepshuffling = false
		for _,pvaldiff := range prevpvaldiffs {
			if pvaldiff > convcond { keepshuffling = true }
		}
		prevpval = pval2
	}
	pval := (pval1 + pval2)/2.0
	return pval
}

func shuffledoms(clussizes []int, r *rand.Rand) ([][]int, []int) {

	permsizes := make([]int, len(clussizes))
	perm := r.Perm(len(clussizes))
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

//	wg.Wait()
	// find how many times shuffvi values are less than given vival
	/*count1 := 0
	count2 := 0
	for i:=0; i< nshuffles; i++ {
		if shuffvi1[i] - vival < 1e-10 { count1++ }
		if shuffvi2[i] - vival < 1e-10 { count2++ }
	}*/
	//pval := (float64(count1+1)/float64(nshuffles+1) + float64(count2+1)/float64(nshuffles+1))/2.0
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
