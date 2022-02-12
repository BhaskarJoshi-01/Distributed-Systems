import sys

ccs = []

for line in sys.stdin:
    cc = set([int(ele) for ele in line.strip().split()])
    for i in range(len(ccs)):
        if len(cc.intersection(ccs[i])) > 0:
            ccs[i] = cc.union(ccs[i])
            break
    else:
        ccs.append(cc)

for cc in ccs:
    for ele in cc:
        print(ele, end=" ")
    print()

        

