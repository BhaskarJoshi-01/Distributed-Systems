#!/usr/bin/python
import os
import sys


def preprocess(init_inp):
    file = open(init_inp, "r")

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

    out_file = open(init_inp, "w+")
    for row_idx, row in enumerate(matA):
        out_file.write(f"A {row_idx} {p} {','.join([str(ele) for ele in row])}\n")

    for row_idx, row in enumerate(matB):
        out_file.write(f"B {row_idx} {m} {','.join([str(ele) for ele in row])}\n")
    
    out_file.close()





if (len (sys.argv) != 7):
    print ("Usage: " + sys.argv[0] + " <jar file address> <input_file> <hdfs_input_dir> <hdfs_output_dir> <mapper_dir> <reducer_dir>")
    exit (1)

jar_file, input_file, hdfs_input_dir, hdfs_output_dir, mapper_dir, reducer_dir = sys.argv[1:]

preprocess(input_file)


# Set uid and gid to the current user
os.setuid(os.getuid())

# Delete both dirs
os.system('hdfs dfs -rm -R ' + hdfs_input_dir)
os.system('hdfs dfs -rm -R ' + hdfs_output_dir)

# create the input dir and copy the file
os.system('hdfs dfs -mkdir ' + hdfs_input_dir)
os.system('hdfs dfs -copyFromLocal ' + input_file + ' ' + hdfs_input_dir)

# get the input file name / or \ or //
input_file =  input_file.split('/')[-1]
input_file = input_file.split('\\')[-1]
input_file = input_file.split('//')[-1]
# Run mapreduce
final_command='hadoop jar ' + jar_file + ' -input ' + hdfs_input_dir + input_file + ' -output ' + hdfs_output_dir + ' -mapper ' + '"python3 ' + mapper_dir + 'mapper.py"' + ' -reducer ' + '"python3 ' + reducer_dir + 'reducer.py"'
print('command: ' + final_command)
os.system(final_command)


# python3 runner.py ../jar_files/hadoop-streaming-3.3.1.jar ./input.txt /Q1/ /output /mnt/Q1/ /mnt/Q1/
