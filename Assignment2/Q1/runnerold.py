file = open("input.txt", "r")

line = [int(ele) for ele in file.readline().strip().split()]
m = line[0]
n = line[1]
matA = []
for i in range(line[0]):
    matA.append([int(ele) for ele in file.readline().strip().split()])

line = [int(ele) for ele in file.readline().strip().split()]
p = line[1]
matB = []
for i in range(line[0]):
    matB.append([int(ele) for ele in file.readline().strip().split()])
file.close()

for row_idx, row in enumerate(matA):
    print("A", row_idx, p, ",".join([str(ele) for ele in row]))

for row_idx, row in enumerate(matB):
    print("B", row_idx, m, ",".join([str(ele) for ele in row]))

