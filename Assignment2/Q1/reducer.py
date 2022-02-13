
#!/usr/bin/python3

# import sys because we need to read and write data to STDIN and STDOUT
import sys

# {MAT}A/B ri p/m vals 
# 
# reading entire line from STDIN (standard input)


ele = 0
i_prev = 0
k_prev = 0
lines = []

for line in sys.stdin:
    lines.append(line)
    if len(lines) < 2:
        continue
    # print(lines)
    line = lines[0].strip().split()
    i, k, pos = (int(ele) for ele in line[0].strip("()").split(","))
    val1 = int(line[1])

    line = lines[1].strip().split()
    i, k, pos = (int(ele) for ele in line[0].strip("()").split(","))
    val2 = int(line[1])

    if (i,k) == (i_prev, k_prev):
        ele += val1 * val2
    else:
        print(ele, end=" ")
        ele = val1 * val2
        if i != i_prev:
            print()
        i_prev = i
        k_prev = k
    lines = []

print(ele)

