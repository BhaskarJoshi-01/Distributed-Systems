import random

edges = int(input("Enter edges: "))

edge_set = set()

for i in range(edges):
    u = random.randint(1, 100)
    v = random.randint(1, 100)
    if(u > v):
        temp = u
        u = v
        v = temp
    edge_set.add((u, v))

file = open("input", "w+")
for edge in edge_set:
    file.write(str(edge[0]) + " " + str(edge[1]) + "\n")
file.close()


#####################################

# Python program to print connected
# components in an undirected graph


class Graph:

	# init function to declare class variables
	def __init__(self, V):
		self.V = V
		self.adj = [[] for i in range(V)]

	def DFSUtil(self, temp, v, visited):

		# Mark the current vertex as visited
		visited[v] = True

		# Store the vertex to list
		temp.append(v)

		# Repeat for all vertices adjacent
		# to this vertex v
		for i in self.adj[v]:
			if visited[i] == False:

				# Update the list
				temp = self.DFSUtil(temp, i, visited)
		return temp

	# method to add an undirected edge
	def addEdge(self, v, w):
		self.adj[v].append(w)
		self.adj[w].append(v)

	# Method to retrieve connected components
	# in an undirected graph
	def connectedComponents(self):
		visited = []
		cc = []
		for i in range(self.V):
			visited.append(False)
		for v in range(self.V):
			if visited[v] == False:
				temp = []
				cc.append(self.DFSUtil(temp, v, visited))
		return cc


# Driver Code

g = Graph(100)
for edge in edge_set:
	g.addEdge(edge[0] - 1, edge[1] - 1)
	
cc_list = g.connectedComponents()
cc_set = []
for cc in cc_list:
    cc_set.append(set([ele +1 for ele in cc]))

cc_set.sort()

####################################3

input("Press Enter to continue...")

of = open("output.txt", "r")
ccs = []
for line in of.readlines():
    cc = set([int(ele) for ele in line.strip().split()])
    ccs.append(cc)

ccs.sort()

def get_set(ccs):
    ret = []
    for cc in ccs:
        st = ""
        for ele in cc:
            st += str(ele) + ","
        ret.append(st)
    ret.sort()
    return ret

print(get_set(ccs) == get_set(cc_set))

