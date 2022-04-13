// SLURM PARAMETERS
// #!/bin/bash
// #SBATCH --error=err-%j.err
// #SBATCH --output=par_task%j.out
// #SBATCH --partition=physical
// #SBATCH --time=1:00:00
// #SBATCH --nodes=1
// #SBATCH --cpus-per-task=32
// #SBATCH --ntasks-per-node=1
// #SBATCH --mem=32G

// module load gcc/10.1.0
// export OMP_PLACES="{0:16:2}, {1:16:2}"
// g++ -fopenmp -Wall -o par_task kseqalign_taskalloc.cpp
// cat mseq.dat | ./par_task

#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>
#include <omp.h>
#include "sha512.hh"
#include <unordered_map>
#include <vector>
#include <bits/stdc++.h>
#include <iterator>
#include <math.h>
#include <queue>

#define MAX_THREAD_NUM 32
#define SYSTEM_MEMORY 29
#define MIN(a,b) (a < b ? a : b)

using namespace std;

string getMinimumPenalties(string *genes, int k, int pxy, int pgap, int *penalties);
uint64_t estimateMemorySize(int i, int j);
int estimateThreadNumber(int i, int j);
int **new2d (int width, int height);
int **generateDPTable(string x, string y, int pxy, int pgap, int numThreads);
int getMinimumPenalty(int **dp, std::string x, std::string y, int pxy, int pgap, int *xans, int *yans);

struct problemInfo {
    string x;
    string y;
    string key;

    uint64_t memorySize;
    int threadNum;
    string hashResult;
    int penalty;
};

void showMap(unordered_map<string, problemInfo> map);
void showVector(vector<pair<int, string>> vect);

uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
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

int main(int argc, char **argv) {
    int misMatchPenalty;
    int gapPenalty;
    int k;

    cin >> misMatchPenalty;
    cin >> gapPenalty;
    cin >> k;

    string genes[k];
    for(int i = 0; i < k; ++i) cin >> genes[i];

    int numPairs = k*(k-1)/2;
    int penalties[numPairs];

    uint64_t start = GetTimeStamp();
    // return all the penalties and the hash of all allignments
    string alignmentHash = getMinimumPenalties(genes, k, misMatchPenalty, gapPenalty, penalties);

    // print the time taken to do the computation
	printf("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));

    // print the alginment hash
	std::cout<<alignmentHash<<std::endl;

	for(int i=0;i<numPairs;i++){
		std::cout<<penalties[i] << " ";
	}
	std::cout << std::endl;
	return 0;
}

string getMinimumPenalties(string *genes, int k, int pxy, int pgap, int *penalties) {
    // Create a map that contains
    unordered_map<string, problemInfo> problemMap;
    vector<pair<uint64_t, string>> memoryVect;

    int completed = 0;
    bool setDynamicTaskAlloc = false;
    int64_t dynamicTaskAllocRequirement = 0.6 * (int64_t)SYSTEM_MEMORY * 1000000000;
    // Populate the map with values
    // Note: the order is essential to ensure that the hash is consistent with sequential solution
    for(int i = 1; i < k; i++) {
        for(int j = 0; j < i; j++) {
            string key = to_string(i) + to_string(j);
            struct problemInfo value;
            value.x = genes[i];
            value.y = genes[j];
            value.key = key;

            // Calculate the memory requirement of the problem
            value.memorySize = estimateMemorySize(value.x.length() + 1, value.y.length() + 1); 

            // To prevent overflows, we check if the requirement is met
            if (!setDynamicTaskAlloc) {
                dynamicTaskAllocRequirement -= value.memorySize;
                if (dynamicTaskAllocRequirement <= 0 || k > 5) {
                    setDynamicTaskAlloc = true;
                }
            } 

            // Calculate the threads required for the problem
            value.threadNum = estimateThreadNumber(value.x.length() + 1, value.y.length() + 1);
            // A place to store the results 
            value.hashResult = "";
            // A place to store the minimum penalty
            value.penalty = 0;

            problemMap[key] = value;
            // Insert the values into the vector
            memoryVect.push_back(make_pair(value.memorySize, value.key));
        }
    }

    // Sort the vector in descending order
    sort(memoryVect.rbegin(), memoryVect.rend());

    // Only do dynamic task allocation when there is a high workload
    if(setDynamicTaskAlloc) {
        // While there are elements in the vector, continuously process the information
        uint64_t availableMemory = SYSTEM_MEMORY;
        availableMemory *= 1000000000;
        int availableThread = MAX_THREAD_NUM;
        int vsize = memoryVect.size();

        while(vsize != completed) {
            string problemKey = "";

            // Check if there are any free threads or free memory to execute task
            if(availableThread > 0 && availableMemory > (uint64_t)memoryVect.back().first && memoryVect.size() > 0) {
                // There is free space and thread to perform calculation
                // Fetch the element that requires the largest memory available
                // We must get a element back because of condition 2
                
                vector<string> execKeyVector = {};

                int vectSize = memoryVect.size();
                int offset = 0;
                for(int i = 0; i < vectSize; ++i) {
                    // Check if the element satisfies the available memory
                    int numElement = i - offset;
                    if(memoryVect.at(numElement).first < availableMemory) {
                        // Check if the element satisfies the thread constraint
                        if(problemMap[memoryVect.at(numElement).second].threadNum + 1 < availableThread) {
                            problemKey = memoryVect.at(numElement).second;
                            availableMemory -= problemMap[problemKey].memorySize;
                            availableThread -= problemMap[problemKey].threadNum + 1;
                            execKeyVector.push_back(problemKey);

                            vector<pair<uint64_t, string>>::iterator it = memoryVect.begin() + i - offset;
                            memoryVect.erase(it);
                            offset += 1;
                        }
                    }
                }

                int taskThread = execKeyVector.size();
                #pragma omp parallel for num_threads(taskThread)
                for(int i = 0; i < taskThread; i++) {
                    printf("%d ", omp_get_num_threads());
                    string problemKey = execKeyVector.at(i);
                    // Create threads to execute stuff here and input results to the map
                    // Generate a table with optimal substructure answers
                    int l = problemMap[problemKey].x.length() + problemMap[problemKey].y.length();
                    int xans[l+1], yans[l+1];
                    int **dp = generateDPTable(problemMap[problemKey].x,
                                            problemMap[problemKey].y,
                                            pxy,
                                            pgap,
                                            problemMap[problemKey].threadNum);

                    #pragma omp atomic
                    availableThread += problemMap[problemKey].threadNum;

                    // Build path from the generated table
                    int penalty = getMinimumPenalty(dp, 
                                                    problemMap[problemKey].x,
                                                    problemMap[problemKey].y,
                                                    pxy,
                                                    pgap,
                                                    xans, 
                                                    yans
                                                    );
                    // Since we have assumed the answer to be n+m long,
                    // we need to remove the extra gaps in the starting
                    // id represents the index from which the arrays
                    // xans, yans are useful
                    int id = 1;
                    int a;
                    for (a = l; a >= 1; a--)
                    {
                        if ((char)yans[a] == '_' && (char)xans[a] == '_')
                        {
                            id = a + 1;
                            break;
                        }
                    }
                    std::string align1="";
                    std::string align2="";
                    for (a = id; a <= l; a++)
                    {
                        align1.append(1,(char)xans[a]);
                    }
                    for (a = id; a <= l; a++)
                    {
                        align2.append(1,(char)yans[a]);
                    }
                    std::string align1hash = sw::sha512::calculate(align1);
                    std::string align2hash = sw::sha512::calculate(align2);
                    std::string problemhash = sw::sha512::calculate(align1hash.append(align2hash));

                    // Uncomment for testing purposes
                    // std::cout << align1 << std::endl;
                    // std::cout << align2 << std::endl;
                    // std::cout << std::endl;

                    // Append the hash into the map
                    problemMap[problemKey].hashResult = problemhash;
                    problemMap[problemKey].penalty = penalty;

                    // Free the memory, thread and remove from vector
                    #pragma omp atomic
                    completed += 1;
                    
                    #pragma omp atomic
                    availableMemory += problemMap[problemKey].memorySize;

                    availableThread += 1;
                }
            }
        }
    } else {
        // In this implementation, no attempts will be made to parallelise the task execution process
        for(int i = 0; i < k; ++i) {
            for(int j = 0; j < i; ++j) {
                string gene1 = genes[i];
                string gene2 = genes[j];

                int m = gene1.length();
                int n = gene2.length();
                int l = m+n;
                int xans[l+1], yans[l+1];

                // Estimate the number of threads needed for the problem
                int threadNum = estimateThreadNumber(m + 1, n + 1);

                int **dp = generateDPTable(gene1, gene2, pxy, pgap, threadNum);
                int penalty = getMinimumPenalty(dp, 
                                                gene1,
                                                gene2,
                                                pxy,
                                                pgap,
                                                xans,
                                                yans);
                // Since we have assumed the answer to be n+m long,
                // we need to remove the extra gaps in the starting
                // id represents the index from which the arrays
                // xans, yans are useful
                int id = 1;
                int a;
                for (a = l; a >= 1; a--)
                {
                    if ((char)yans[a] == '_' && (char)xans[a] == '_')
                    {
                        id = a + 1;
                        break;
                    }
                }
                std::string align1="";
                std::string align2="";
                for (a = id; a <= l; a++)
                {
                    align1.append(1,(char)xans[a]);
                }
                for (a = id; a <= l; a++)
                {
                    align2.append(1,(char)yans[a]);
                }
                std::string align1hash = sw::sha512::calculate(align1);
                std::string align2hash = sw::sha512::calculate(align2);
                std::string problemhash = sw::sha512::calculate(align1hash.append(align2hash));

                // Uncomment for testing purposes
                // std::cout << align1 << std::endl;
                // std::cout << align2 << std::endl;
                // std::cout << std::endl;

                // Append the hash into the map
                string key = to_string(i) + to_string(j);
                problemMap[key].hashResult = problemhash;
                problemMap[key].penalty = penalty;
            }
        }
    }

    // compute the alignment hash from the unordered map
    string alignmentHash = "";
    int probNum = 0;
    for(int i = 1; i < k; ++i) {
        for(int j = 0; j < i; ++j) {
            string key = to_string(i) + to_string(j);

            // Get the problemHash from the map
            alignmentHash = sw::sha512::calculate(alignmentHash.append(problemMap[key].hashResult));

            // Update penalties
            penalties[probNum] = problemMap[key].penalty;
            probNum++;
        }
    }

    // release the memory held by the map and vectors
    problemMap.clear();
    memoryVect.clear();

    return alignmentHash;
}

// Function to display all elements within the problem map
// Note: This is for DEBUGING only
void showMap(unordered_map<string, problemInfo> map) {
    for (auto it : map) {
        cout << "Information about: " << it.second.key << endl;
        cout << "String 1: " << it.second.x << endl;
        cout << "String 2: " << it.second.y << endl;
        cout << "MemorySize: " << it.second.memorySize << endl;
        cout << "threadNum: " << it.second.threadNum << endl;
        cout << "" << endl;
    }
}

// Function to display all the elements within a vector
// Note: This is for debugging only
void showVector(vector<pair<int, string>> vect) {
    for (auto i : vect) {
        cout << i.second << endl;
    }
}

uint64_t estimateMemorySize(int i, int j) {
    // Assume that the table is in the form of int
    return i * j * sizeof(int);
}

int estimateThreadNumber(int i, int j) {
    // We have problem size 10000 vs 2000 tile size as the baseline
    // Experimentally, we have determined that 3/2 is an ideal constant for 
    // the ratio of increase/decrese of computing power against the increase/
    // decrese of tile size
    // First, we determine the r value
    float r1 = (i > 10000) ? 2 : 0.5;
    float r2 = (i > 10000) ? 1.5 : 2/3;

    // calculate the tile size 
    int tileI;
    tileI = (int)2000*r2*(float)(log(( (float) i/10000))/log(r1));

    // to handle cases where tileI is 0
    tileI = (tileI == 0) ? 2000 : tileI;

    int threadNumI = i/tileI;
    // Calculate the r values for j
    r1 = (j > 10000) ? 2 : 0.5;
    r2 = (j > 10000) ? 1.5 : 2/3;

    // calculate the tile size
    int tileJ;
    
    tileJ = (int)2000*r2*(float)(log((float) (j/10000))/log(r1));
    
    // to handle cases where tileI is 0
    tileJ = (tileJ == 0) ? 2000 : tileJ;
    
    int threadNumJ = j / tileJ;
    int threadNum =  MIN(threadNumI, threadNumJ);
    return (threadNum == 0) ? threadNum + 1 : threadNum;
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

// Generate the optimal substructure table
int **generateDPTable(string x, string y, int pxy, int pgap, int numThreads) {
    int m = x.length();
    int n = y.length();

    int di = (m+1)/numThreads;
    int dj = (n+1)/numThreads;

    int i = 0;
    int j = 0;

    int height = m + 1;
    int width = n + 1;

    // Create a 2d table for storing the optimal substructure answers
    int **dp = new2d(height, width);

    int diagonal = (height/di) + (width/dj);
    diagonal += ((height % di) + (width % dj) > 0) ? 1 : 0;

    // Separate parallelised code section and sequential code section
    // There will be sections where sequential code is used as the problem size
    // may be too small to justify the use of parallelism
    if(numThreads > 1) {
        // Parallelised table creation to utilise the first touch memory policy
        for(int d = 0; d < diagonal; ++d) {
            int diff_i, diff_j, diag_i, diag_j, length, tile;

            diff_i = height - i - 1;
            diff_j = j;

            diag_i = 1 + ((diff_i) / di);
            diag_j = 1 + ((diff_j) / dj);
            length = MIN(diag_i, diag_j);

            // We parallelise the tile initialisation to set the memory so that threads
            // close to it executes it.
            #pragma omp parallel for private(tile) schedule(static) num_threads(numThreads)
            for(tile = 0; tile < length; ++tile) {
                int ii = i + (tile * di);
                int jj = j - (tile * dj);

                int imax = MIN(ii + di, height);
                int jmax = MIN(jj + dj, width);

                for(int iii = ii; iii < imax; ++iii) {
                    for(int jjj = jj; jjj < jmax; ++jjj) {
                        dp[iii][jjj] = omp_get_thread_num();
                    }
                }
            }

            j += dj;
            if ( j >= width) {
                j -= dj;
                i += di;
            }
        }

        // Initialise the table 
        for(i = 0; i <= m; ++i) {
            dp[i][0] = i*pgap;
        }

        for(i = 0; i <= n; ++i) {
            dp[0][i] = i*pgap;
        }

        // reset the value of i and j
        i = 0, j = 0;

        // Apply the substructure optimal answers onto the table
        for(int d = 0; d < diagonal; ++d) {
            int diff_i, diff_j, diag_i, diag_j, length, tile;

            diff_i = height-i-1;
            diff_j = j;

            diag_i = 1 + ((diff_i) / di);
            diag_j = 1 + ((diff_j) / dj);
            length = MIN(diag_i, diag_j);

            // We parallelise the tile calculation
            #pragma omp parallel for private(tile) schedule(static) num_threads(numThreads) proc_bind(close)
            for(tile = 0; tile < length; ++tile) {
                int ii = i + (tile * di);
                int jj = j - (tile * dj);

                int imax = MIN(ii + di, height);
                int jmax = MIN(jj + dj, width);

                for(int iii = ii; iii < imax; ++iii) {
                    for(int jjj = jj; jjj < jmax; ++jjj) {
                        // Since the edges of the table has been initialised, no
                        //calculations needed
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
            j += dj;
            if (j >= width) {
                j -= dj;
                i += di;
            }
        }
    } else {
        // Sequential code
        // Set the values for the table
        size_t size = m + 1;
        size *= n + 1;

        // intialising the table
        for (i = 0; i <= m; i++)
        {
            dp[i][0] = i * pgap;
        }
        for (i = 0; i <= n; i++)
        {
            dp[0][i] = i * pgap;
        }

        // calcuting the minimum penalty
        for (i = 1; i <= m; i++)
        {
            for (j = 1; j <= n; j++)
            {
                if (x[i - 1] == y[j - 1])
                {
                    dp[i][j] = dp[i - 1][j - 1];
                }
                else
                {
                    dp[i][j] = min3(dp[i - 1][j - 1] + pxy ,
                            dp[i - 1][j] + pgap ,
                            dp[i][j - 1] + pgap);
                }
            }
        }
    }

    return dp;
}

// Function to find the alignment of the two sequences given that the optimal substructure table is generated
int getMinimumPenalty(int **dp, std::string x, std::string y, int pxy, int pgap, int *xans, int *yans) {
    int m = x.length();
    int n = y.length();

    // Reconstruct the solution
    int l = n + m; // maximum possible length
	
	int i = m; 
    int j = n;
	
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