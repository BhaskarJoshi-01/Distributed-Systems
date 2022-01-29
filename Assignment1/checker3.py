import os
import sys
import numpy as np

# with open(sys.argv[1], 'r') as f:
#     n = int(f.readline())
#     data = f.readlines()
#     data = [[float(word) for word in line.split()] for line in data][:n]
#     data = np.array(data)
#     # print(data)
#     A = data[:n, :n]
#     B = data[:n, n]

def gen_random_data(size):
    while 1:
        a = np.random.randint(-1e7, 1e7, size*size).reshape(size, size)
        if np.linalg.norm(a) > 0:
            return a, np.random.randint(-1e7, 1e7, size)


def get_my_ans():
    os.system("mpirun -np 11  --oversubscribe ./a.out in.txt o3_1.txt") #type the command that u use to compile
    with open("o3_1.txt") as f: #type the name of the output file
        data = f.readlines()[0]
    return [float(i) for i in data.split()]

for i in np.random.randint(1, 101, 1000):
    A,B = gen_random_data(i)
    ans = ((np.linalg.inv(A)@B))
    with open('in.txt', 'w') as f:
        s = str(i)+"\n"
        tmp = np.c_[A, B]
        for row in tmp:
            for val in row:
                s+= str(val) + " "
            s+= "\n"
        
        f.write(s)
    my_ans = get_my_ans()
    print(my_ans, ans)
    assert(np.isclose(get_my_ans(), ans).all())



