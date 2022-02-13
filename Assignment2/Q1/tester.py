import numpy as np

m = int(input("m: "))
n = int(input("n: "))
p = int(input("p: "))

matA = np.random.randint(-1000, 1000, (m, n))
matB = np.random.randint(-1000, 1000, (n, p))

ifile = open("input.txt", "w+")
ifile.write(f"{m} {n}\n")
for row in matA:
    ifile.write(f"{' '.join([str(ele) for ele in row])}\n")
ifile.write(f"{n} {p}\n")
for row in matB:
    ifile.write(f"{' '.join([str(ele) for ele in row])}\n")
ifile.close()


res = np.matmul(matA, matB)

input("Press Enter to continue...")

ofile = open("output.txt", "r")

resmat = []
for i in range(m):
    resmat.append([int(ele) for ele in ofile.readline().strip().split()])
ofile.close()

resmat = np.array(resmat)
print((resmat == res).all())
