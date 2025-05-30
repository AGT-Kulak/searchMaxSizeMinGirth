import sys
import os
from subprocess import run, Popen, PIPE
from common import *
import time


if __name__ == "__main__":
   
    iter = int(sys.argv[1])
    L = int(sys.argv[2])
    max_processes = int(sys.argv[3])
    
    file_dir = f"downRun_iter{iter}_{L}"

    bestGraphsSaveDir = f"bestGraphs/iter{iter}Down_{L}"

    try:
        os.makedirs(file_dir)
        os.makedirs(file_dir + "/ptg")
        os.makedirs(file_dir + "/ls")
        os.makedirs(bestGraphsSaveDir)
    except FileExistsError:
        # directory already exists
        pass
 
    readStartGraphs(f"bestGraphs/iter{iter}Up_{L}")

    lowerOrder = 50
    upperOrder = 59

    orders = range(upperOrder-1,lowerOrder-1,-1)
    numItersFunc = lambda x : x * 1000
    p = 0.5
    numDelEdge = 3
    numsDelEdgeFuncs = [lambda x : 3, lambda x : max([x//10,3])]
    input_file_ptg = startGraphFileOfOrder(upperOrder)

    # for local search: take last g6 string in stderr output

    for n in orders:

        bitsetSize = getBitsetSize(n)

        numIters = numItersFunc(n)

        numsDelEdges = set()
        for numsDelEdgeFunc in numsDelEdgeFuncs:
            numsDelEdges.add(numsDelEdgeFunc(n))

        output_file_ptg = f"{file_dir}/ptg/n{n}.g6"

        command = (f"./lib/nauty/delptg -n1 < {input_file_ptg} | ./lib/topLMostEdges {L} > {output_file_ptg}")

        t1 = time.time()
        data = run(command, capture_output=True, shell=True, text=True)
        t2 = time.time()
        print(data.stderr)
        print(f"ptg time = {t2-t1}s")

        print(f"Delete vertex n={n+1}->{n} done")

        numGraphsOutPtg = getNumLinesInFile(output_file_ptg)

        processes = []  # List to store process objects
        curr_num_processes = 0
        finished = [False] * numGraphsOutPtg * 2

        output_files_ls = []
        for lineNum in range(1,numGraphsOutPtg+1):

            for numDelEdge in numsDelEdges:

                recency=-1
                if numDelEdge == 3:
                    recency = n
                else:
                    recency = n//3

                output_file_ls = f"{file_dir}/ls/out_ls_n{n}_lineNum{lineNum}_p{p}_recency{recency}_numDel{numDelEdge}_iters{numIters}.g6"
                command=(f"sed '{lineNum}!d' {output_file_ptg} | ./searchMaxSizeMinGirth{bitsetSize} -n{n} -r -p{p} -R{recency} -g5 -d{numDelEdge} -i{numIters} 1>{output_file_ls}")
                print(command)

                # Wait if we have reached the maximum allowed concurrent processes.
                while curr_num_processes >= max_processes:
                    # Check for any process that has completed.
                    for i,process in enumerate(processes):
                        if process.poll() is not None and not finished[i]:  # Process finished
                            curr_num_processes -= 1
                            finished[i] = True
                        
                    # If none have finished, wait a short time before checking again.
                    time.sleep(2)

                # Start the process in parallel, capturing stderr.
                process = Popen(command, shell=True, text=True, stderr=PIPE)
                processes.append(process)
                output_files_ls.append(output_file_ls)
                curr_num_processes += 1

        for process in processes:
            _, stderr = process.communicate()  # Wait for process to complete and get stderr
            if stderr:
                print(stderr)  # Print stderr if there's any error output


        print(f"ls n={n} done")
    

        if isStartGraphOfOrder(n):
            output_files_ls.append(startGraphFileOfOrder(n))

        input_file_ptg = f"{bestGraphsSaveDir}/n{n}.g6"

        shortg_command = "cat "
        for output_file_ls in output_files_ls:
            shortg_command += output_file_ls + " "
        shortg_command += f"| ./lib/nauty/shortg | ./lib/topLMostEdges {L} > {input_file_ptg}"
        print(shortg_command)
        t1 = time.time()
        data = run(shortg_command, capture_output=True, shell=True, text=True)
        t2 = time.time()
        print(data.stderr)
        print(f"cat time = {t2-t1}s")

        maxEdgesLs = getMostEdges([input_file_ptg])

        print()
        print(f"Best num edges n={n}: {maxEdgesLs}")
        print()

        


    
