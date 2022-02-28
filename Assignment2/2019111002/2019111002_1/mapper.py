
#!/usr/bin/python3

# import sys because we need to read and write data to STDIN and STDOUT
import sys

# {MAT}A/B ri p/m vals 
# 
# reading entire line from STDIN (standard input)

def str_num(num):
    ret = str(num)
    while len(ret) != 3:
        ret = "0" + ret
    return ret

for line in sys.stdin:
    line = line.strip().split()
    mat = line[0]
    ri = int(line[1])
    p_or_m = int(line[2])
    vals = [int(ele) for ele in line[3].strip().split(',')]
    if mat == "A":
        for j, ele in enumerate(vals):
            for k in range(p_or_m):
                print(f"({str_num(ri)},{str_num(k)},{str_num(j)}) {ele}")
    else:
        for k, ele in enumerate(vals):
            for i in range(p_or_m):
                print(f"({str_num(i)},{str_num(k)},{str_num(ri)}) {ele}")
