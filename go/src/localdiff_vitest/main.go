package main

import (
	"math/rand"
	"strconv"
	"strings"
	"fmt"
	"hicutil"
	"math"
//	"os"
)



func main() {

	filelist := []string{"/mnt/disk17/armatusresults/A549/100kb/EncodeChr_A549_combo_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/Caki2/100kb/EncodeChr_Caki2_combo_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/G401/100kb/EncodeChr_G401_combo_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/LNCaP-FGC/100kb/EncodeChr_LNCaP-FGC_combo_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/NCI-H460/100kb/EncodeChr_NCI-H460_combo_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/Panc1/100kb/EncodeChr_Panc1_combo_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/RPMI-7951/100kb/EncodeChr_RPMI-7951_combo_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/SJCRH30/100kb/EncodeChr_SJCRH30_combo_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/SKMEL5/100kb/EncodeChr_SKMEL5_combo_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/SKNDZ/100kb/EncodeChr_SKNDZ_combo_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/SKNMC/100kb/EncodeChr_SKNMC_combo_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/T47D/100kb/EncodeChr_T47D_combo_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/IMR90/100kb/RaoChr_IMR90_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/GM12878/100kb/RaoChr_GM12878_rep1_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/HMEC/100kb/RaoChr_HMEC_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/HUVEC/100kb/RaoChr_HUVEC_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/K562/100kb/RaoChr_K562_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/KBM7/100kb/RaoChr_KBM7_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/NHEK/100kb/RaoChr_NHEK_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/IMR90-Dixon/hicpro/100kb/DixonChr_IMR90_combined_100kb_gmax1.0", "/mnt/disk17/armatusresults/hESC/hicpro/100kb/DixonChr_hESC_combined_100kb_gmax1.0", "/mnt/disk17/armatusresults/GM06990/100kb/LiebAidChr_GM06990_IC_100kb_gmax1.0", "/mnt/disk17/armatusresults/K562/100kb/LiebAidChr_K562_IC_100kb_gmax1.0"}

	res := 100000
	for i:=0; i < 10; i++{
		celltype1 := rand.Intn(len(filelist))
		celltype2 := rand.Intn(len(filelist))
		for celltype1 == celltype2 {
			celltype2 = rand.Intn(len(filelist))
		}
		
		chrnum := strconv.Itoa(rand.Intn(21)+1)
		
		filename1 := filelist[celltype1]
		idx := strings.Index(filename1, "Chr")
		filename1 = filename1[:idx+3]+chrnum+filename1[idx+3:]

		filename2 := filelist[celltype2]
		idx = strings.Index(filename2, "Chr")
		filename2 = filename2[:idx+3]+chrnum+filename2[idx+3:]
		
		//read in TAD files
		tadlists := processTADLists([]string{filename1,filename2}, &res, true, 8.8)
		//fmt.Println(tadlists)

		// calculate VI matrix
		var sizeM int
		if tadlists[0][len(tadlists[0])-1][1] > tadlists[1][len(tadlists[1])-1][1] {
			sizeM = tadlists[0][len(tadlists[0])-1][1]
		} else {
			sizeM = tadlists[1][len(tadlists[1])-1][1]
		}
		M := make([][]float64, sizeM)//should be max(tadlists)
		for i := range M {
			M[i] = make([]float64, sizeM)
			for j := i+1; j < sizeM; j++ {
				// calculate VI
				n := j-i+1
				intvl1 := hicutil.ProcessIntervals(tadlists[0], i, j)
				intvl2 := hicutil.ProcessIntervals(tadlists[1], i, j)
				if (intvl1[0][0] == 0 && intvl1[0][1] == n-1) || (intvl2[0][0] == 0 && intvl2[0][1] == n-1) {
					M[i][j] = math.NaN()
				} else {
					overlaps := hicutil.CalcOverlaps(intvl1, intvl2)
					clus1sizes := make([]int, len(intvl1))
					for c,clus := range intvl1 {
						clus1sizes[c] = clus[1]-clus[0]+1 }
					clus2sizes := make([]int, len(intvl2))
					for c,clus := range intvl2 {
						clus2sizes[c] = clus[1]-clus[0]+1 }
					condh1 := hicutil.CalcCondEntropy(transpose(overlaps), clus1sizes, n)
					condh2 := hicutil.CalcCondEntropy(overlaps, clus2sizes, n)
					M[i][j] = (condh1 + condh2)/math.Log(float64(n))
				}
				//M[i][j] = (condh1 + condh2)/math.Log(float64(n))
//				fmt.Println(M[i][j])
			}
		}
		
		//find local mins
		sizeM += -1
		var localmins [][]int
		atbdy := 0
		for i,row := range M {
			for j,val := range row {
				if j < i+1 { continue }
				// otherwise check if val < all 8 neighbors
				min := true
				if math.IsNaN(val) { 
					min = false
				}
				var neighbors []float64
				if i > 0 && i < sizeM && j < sizeM { //not on an edge
					neighbors = []float64{M[i-1][j-1], M[i-1][j], M[i-1][j+1], M[i][j-1], M[i][j+1], M[i+1][j-1], M[i+1][j], M[i+1][j+1]}
				} else if j == sizeM && i > 0 && i < sizeM {
					neighbors = []float64{M[i-1][j-1], M[i-1][j], M[i][j-1], M[i+1][j-1], M[i+1][j]}
				} else if i == 0 && j < sizeM {
					neighbors = []float64{M[i][j-1],M[i][j+1],M[i+1][j-1], M[i+1][j], M[i+1][j+1]}
				} else if i == sizeM && j < sizeM {
					neighbors = []float64{M[i-1][j-1], M[i-1][j], M[i-1][j+1], M[i][j-1], M[i][j+1]}
				} else if i == 0 && j == sizeM {
					neighbors = []float64{M[i][j-1], M[i+1][j-1], M[i+1][j]}
				}
				for _,compval := range neighbors {
					if compval <= val {
						min = false
						break
					}
				}
				if min {
					localmins = append(localmins, []int{i,j})
					// check if both i and j are at a boundary in either TAD set
					if (contains(tadlists[0], i) || contains(tadlists[1],i)) && (contains(tadlists[0],j) || contains(tadlists[1], j)) {
						atbdy +=1
					}
				}
				
			}
		}
		fmt.Println(filename1, filename2)
		//fmt.Println("localmins =", localmins)
		//fmt.Println("tadlist 1 =",tadlists[0])
		//fmt.Println("tadlist 2 =", tadlists[1])
		fmt.Println("Percentage of local mins at boundary points =",float64(atbdy)/float64(len(localmins)))
	}
}


func processTADLists(tadfilelist []string, res *int, gammaopt bool, medtadlen float64) ([][][]int) {
 
         tadlists := make([][][]int, 2)
         chrlength := 0
         var gamma float64
         for i:=0; i < 2; i++ {
                 if gammaopt == true {
                         tadlists[i],gamma = hicutil.ChooseGamma(medtadlen, tadfilelist[i], *res)
                         _ = gamma
                 } else {
                         tadlists[i] = hicutil.ReadTADFile(tadfilelist[i], *res)
                 }
                 n := tadlists[i][len(tadlists[i])-1][1]
                 if chrlength < n+1 {
                         chrlength = n+1
                 }
         }
         for i:=0; i < 2; i++ {
                 tadlists[i] = hicutil.FillinTADList(tadlists[i], chrlength)
         }
 
         return tadlists
 
 }


func transpose(a [][]int) [][]int {
        n := len(a)
        m := len(a[0])
        b := make([][]int, m)
        for j := 0; j < m; j++ {
                b[j] = make([]int, n)
        }
        for i := 0; i < n; i++ {
                for j := 0; j < m; j++ {
                        b[j][i] = a[i][j]
                }
        }
        return b
}

func contains(a [][]int, x int) bool {
        for _, row := range a {
                if row[1] == x || row[0] == x {
                        return true
                }
        }
        return false
}
