#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>
#include <unordered_map>
#include "sha512.hh"
#include "verbose.hh"
#include <stdbool.h>

#define MIN(a, b) (a < b ? a : b)
#define MAX(a, b) (a > b ? a : b)

using namespace std;

std::string getMinimumPenalties(std::string *genes, int k, int pxy, int pgap, int *penalties);
void showTable(int **dp, int m, int n);
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans);
void tileProcessingSeq(int **dp, int i, int j, int di, int dj, int tile, int height, int width, std::string x, std::string y, int pxy, int pgap);
void tileProcessingParallel(int **dp, int i, int j, int di, int dj, int tile, int height, int width, std::string x, std::string y, int pxy, int pgap);

uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

// equivalent of  int *dp[width] = new int[height][width]
// but works for width not known at compile time.
// (Delete structure by  delete[] dp[0]; delete[] dp;)
int **new2d (int width, int height)
{
	int **dp = new int *[width];
	size_t size = width;
	size *= height;
	int *dp0 = new int [size];
	if (!dp || !dp0)
	{
	    std::cerr << "getMinimumPenalty: new failed" << std::endl;
	    exit(1);
	}
	dp[0] = dp0;
	for (int i = 1; i < width; i++)
	    dp[i] = dp[i-1] + height;

	return dp;
}

int main(int argc, char **argv) {
    int misMatchPenalty;
	int gapPenalty;
	int k;
	std::cin >> misMatchPenalty;
	std::cin >> gapPenalty;
	std::cin >> k;	
	std::string genes[k];
	for(int i=0;i<k;i++) std::cin >> genes[i];

    int numPairs= k*(k-1)/2;

	int penalties[numPairs];

    uint64_t start = GetTimeStamp ();

    std::string alignmentHash = getMinimumPenalties(genes, k, misMatchPenalty, gapPenalty, penalties);

    printf("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));

    std::cout<<alignmentHash<<std::endl;

    for(int i=0;i<numPairs;i++){
		std::cout<<penalties[i] << " ";
	}
    
	std::cout << std::endl;
	return 0;
}

int min3(int a, int b, int c) {
	if (a <= b && a <= c) {
		return a;
	} else if (b <= a && b <= c) {
		return b;
	} else {
		return c;
	}
}

std::string getMinimumPenalties(std::string *genes, int k, int pxy, int pgap, int *penalties) {
    // int probNum = 0;
    std::string alignmentHash="";

    // Create a map for the alignment hashes to be mapped
    unordered_map<string, string> aHashMap;

    for (int i = 1; i < k; i++) {
        for (int j = 0; j < i; j++) {
            std::string gene1 = genes[i];
            std::string gene2 = genes[j];

            int m = gene1.length();
            int n = gene2.length();
            int l = m + n;
            int xans[l+1], yans[l+1];
            // Generate a table of penalties that can be incurred
            getMinimumPenalty(gene1, gene2, pxy, pgap, xans, yans);

            int id = 1;
            int a;
            for(a = l; a >= 1; a--)
            {
                if ((char)yans[a] == '_' && (char)xans[a] == '_')
                {
                    id = a + 1;
                    break;
                }
            }

            std::string align1="";
            std::string align2="";
            
            for(a = id; a <= l; a++)
            {
                align1.append(1, (char)xans[a]);
            }
            
            for(a = id; a <= l; a++)
            {
                align2.append(1, (char)yans[a]);
            }

            std::string align1hash = sw::sha512::calculate(align1);
            std::string align2hash = sw::sha512::calculate(align2);
            std::string problemhash = sw::sha512::calculate(align1hash.append(align2hash));
            std::string hashKey = align1+align2;

            // TODO: I can only have 1 thread writing at a time
            #pragma omp critical
            {
                aHashMap[hashKey] = problemhash;
            }

            // Uncomment for testing purposes
            // std::cout << penalties[probNum] << std::endl;
            // std::cout << align1 << std::endl;
            // std::cout << align2 << std::endl;
            // std::cout << std::endl;

			// probNum++;
        }
    }

    // Print out all the values in umap
    // unordered_map<string, double>:: iterator itr;
    for(auto x : aHashMap) 
    {
        std::cout << x.first << " " << x.second << endl;
    }

    return alignmentHash;
}

int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans) {
    int diagonal, i, j;

    int m = x.length();
    int n = y.length();

    // Create a 2d table for storing the optimal substructure answers
    int **dp = new2d (m+1, n+1);
	size_t size = m + 1;
	size *= n + 1;
	memset (dp[0], 0, size);

    // Initialising the table
    for(i = 0; i <= m; ++i) {
        dp[i][0] = i * pgap;
    }

    for(i = 0; i <= n; ++i) {
        dp[0][i] = i * pgap;
    }

    int width = n+1;
    int height = m+1;

    int di = 10;
    int dj = 10;

    i = 0, j = 0;

    diagonal = (height/di) + (width/dj);
    diagonal += ((height % di) + (width % dj) > 0) ? 1 : 0;
    
    // #pragma omp parallel for 
    for(int d = 0; d < diagonal; ++d) {
        int diff_i, diff_j, diag_i, diag_j, length, tile;

        diff_i = height-i-1;
        diff_j = j;

        diag_i = 1 + ((diff_i) / di);
        diag_j = 1 + ((diff_j) / dj);
        length = MIN(diag_i, diag_j);

        // We parallelise the tile calculation
        // #pragma omp parallel for
        for(tile = 0; tile < length; ++tile) {
            tileProcessingParallel(dp, i, j, di, dj, tile, height, width, x, y, pxy, pgap);
        }
        j += dj;
        if (j >= width) {
            j -= dj;
            i += di;
        }
    }

    showTable(dp, m+1, n+1);

    // Reconstructing the solution
	int l = n + m; // maximum possible length
	
	i = m; j = n;
	
	int xpos = l;
	int ypos = l;
	
	while ( !(i == 0 || j == 0))
	{
		if (x[i - 1] == y[j - 1])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[i - 1][j - 1] + pxy == dp[i][j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[i - 1][j] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)'_';
			i--;
		}
		else if (dp[i][j - 1] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)'_';
			yans[ypos--] = (int)y[j - 1];
			j--;
		}
	}
	while (xpos > 0)
	{
		if (i > 0) xans[xpos--] = (int)x[--i];
		else xans[xpos--] = (int)'_';
	}
	while (ypos > 0)
	{
		if (j > 0) yans[ypos--] = (int)y[--j];
		else yans[ypos--] = (int)'_';
	}

	int ret = dp[m][n];

	delete[] dp[0];
	delete[] dp;

    return ret;
}

// This function processes each tile sequentially
void tileProcessingSeq(int **dp, int i, int j, int di, int dj, int tile,
                        int height, int width, std::string x, std::string y, int pxy, int pgap) {
    int ii = i + (tile * di);
    int jj = j + (tile * dj);

    int imax = MIN(ii + di, height);
    int jmax = MIN(jj + dj, width);

    for(int iii = ii; iii < imax; ++iii) {
        for(int jjj = jj; jjj < jmax; ++jjj) {
            // Ensure that the edge values are not overwritten
            if(iii == 0 || jjj == 0) {
                // do nothing
            } else {
                if(x[iii - 1] == y[jjj - 1]) {
                    dp[iii][jjj] = dp[iii-1][jjj-1];
                } else {
                    if(x[iii-1] == y[jjj-1]) {
                        dp[iii][jjj] = dp[iii-1][jjj-1];
                    } else {
                        dp[iii][jjj] = min3(dp[iii-1][jjj-1] + pxy,
                                            dp[iii-1][jjj] + pgap,
                                            dp[iii][jjj-1] + pgap);
                    }
                }
            }
        }
    }
}

void tileProcessingParallel(int **dp, int i, int j, int di, int dj, int tile,
                            int height, int width, std::string x, std::string y, int pxy, int pgap) {
    int ii = i + (tile * di);
    int jj = j - (tile * dj);

    int imax = MIN(ii + di, height);
    int jmax = MIN(jj + dj, width);

    // First, get the start point of every iteration
    int tileIterations = MIN(MIN(di, imax - ii), MIN(dj, jmax - jj));
    int iter;

    for(iter = 0; iter < tileIterations; iter++) {
        int iii, jjj;
        int startI = ii + iter;
        int startJ = jj + iter;

        // Check bounds
        if(startI >= imax) {
            startI = imax - 1;
        }

        if(startJ >= jmax) {
            startJ = jmax - 1;
        }

        // Expand the loop length and height wise
        for(iii = startI; iii < imax; iii++) {
            if(iii == 0 || startJ == 0) {
                // do nothing
            } else {
                if(x[iii - 1] == y[startJ - 1]) {
                    dp[iii][startJ] = dp[iii-1][startJ-1];
                } else {
                    dp[iii][startJ] = min3(dp[iii-1][startJ-1] + pxy,
                                        dp[iii-1][startJ] + pgap,
                                        dp[iii][startJ-1] + pgap);
                }
            }
        }

        for(jjj = startJ; jjj < jmax; jjj++) {
            if(jjj == 0 || startI == 0) {
                // do nothing
            } else {
                if(x[startI - 1] == y[jjj - 1]) {
                    dp[startI][jjj] = dp[startI - 1][jjj - 1]; 
                } else {
                    dp[startI][jjj] = min3(dp[startI-1][jjj-1] + pxy,
                                           dp[startI-1][jjj] + pgap,
                                           dp[startI][jjj-1] + pgap);
                }
            }
        }
    }
}

void showTable(int **dp, int m, int n) {
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            printf("%3d", dp[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}