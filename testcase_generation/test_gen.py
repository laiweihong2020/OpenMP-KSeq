import sys
from random import random

def main(args):
    misMatchPenalty = args[1]
    gapPenalty = args[2]
    k = args[3]
    length = args[4]
    filename = args[5]

    outputStr = list()

    # Gene string generation
    for x in range(int(k)):
        str = ""
        for len in range(int(length)):
            value = random()*4
            if(value < 1):
                str += "A"
            elif(value < 2 and value >= 1):
                str += "C"
            elif(value < 3 and value >= 2):
                str += "G"
            else:
                str += "T"
        outputStr.append(str)

    # Write everything to file
    f = open("../testcase/" + filename, "a")
    f.write(misMatchPenalty + "\n")
    f.write(gapPenalty + "\n")
    f.write(k + "\n")
    for x in outputStr:
        f.write(x + "\n")

if __name__ == "__main__":
    args = sys.argv
    main(args)