
#!/usr/bin/python3

# import sys because we need to read and write data to STDIN and STDOUT
import sys
from collections import defaultdict

# reading entire line from STDIN (standard input)
graph = defaultdict(lambda: []) # key: node value: adj list

for line in sys.stdin:

    tokens = line.strip().split()
    u = tokens[0]
    v = tokens[1]
    graph[u].append(v)
    graph[v].append(u)

parent = {}
for node in graph.keys():
    parent[node] = -1

def dfs(node, tag, parent, graph):
    if parent[node] != -1:
        return
    parent[node] = tag
    for child in graph[node]:
        dfs(child, tag, parent, graph)

ccc = 0
for node in graph.keys():
    if parent[node] == -1:
        dfs(node, ccc, parent, graph)
        ccc += 1

for tag in range(ccc):
    for node in graph.keys():
        if parent[node] == tag:
            print(node, end=" ")
    print()





