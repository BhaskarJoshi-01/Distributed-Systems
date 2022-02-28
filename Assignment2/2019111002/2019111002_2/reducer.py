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

all_nodes = set()

for cc in ccs:
    for ele in cc:
        all_nodes.add(ele)
        print(ele, end=" ")
    print()

for node in range(1, 101):
    if node not in all_nodes:
        print(node)

        

