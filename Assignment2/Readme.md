# Assignment 2 
# Map Reduce: Distributed Systems
## Bhaskar Joshi
## 2019111002

# Describing how the environment was setup 
I have attempted the assignment in Hadoop that was basically the requirement and it was set up in Docker. The important aspect of this was I was able to make updates in one file and put that in shared component and it was barely a script's work to make changes and it was working perfectly fine.
The instructions to setup Docker can be found on hdoop folder `https://github.com/BhaskarJoshi-01/Distributed-Systems/tree/main/Assignment2/hdoop `
Corresponding jar files can be found on `https://github.com/BhaskarJoshi-01/Distributed-Systems/tree/main/Assignment2/jar_files`

## Q1
  
- Problem statement:  Given two matrices A of size m ∗ n and B of size n ∗ p. Output the matrix
multiplication of A and B
- Input: The first line of the input contains m and n followed by m lines of elements
belonging to matrix A. Then it is followed by n and p, followed by n lines of
elements belonging to matrix B.
- Output: The output should contain lines of elements belonging to output A ∗ B.

Runner Script:
This file reads input from `input.txt file`, hopefully given in the above format.
What it does is splits the input into m,n,p, A and B mentioned as per question above.
Then from the inputs I am sending :
```Matrix row_index (k/i) values ```
Let me give an example , we have A= [[1,2],[3,4]] B=[[2,4],[4,8]]
then we will have 
`Matrix(A/B) row index (p/m) row values ` 

```
A 0 2 1,2
A 1 2 3,4
B 0 2 2,4
B 1 2 4,8
```
And this was `redirected to input file` (can change later) for mapper.

Mapper.py:
After we get our values in `Matrix(A/B) row index (p/m) row values ` format , we try to format it to `(Row_index,P,j) , value of corresponding row ` for Matrix A and `(M,value,Row_index) , value of corresponding row ` for Matrix B.
Basically the first term in the bracket above is the key and it is basically (i,k,j) and its corresponding value for matrix A like the values assoicated and its positions.
The problem that I faced was since it was a key and so has to be sorted but Ascii value of `,` was less than numbers so it was not getting properly ordered.
The solution I came up with was append 0 in the begining of the number if there were less than 3 digits in number making use of the constraints provided for n,m,p values. And so the data was processed on mapper and was waiting to be sent to reducer for further processing.

Reducer.py:
The key,value obtained from mapper was taken and what I did was:
for every j :
take submission of every possible A_ij * B_jk 
and what we recieve from this is for every (i,k) the corresponding output value.

I have also made tester script for checking A,B upto 999*999.Cause it was minimun requirement and this could be increased/decreased as per requirement from  https://github.com/BhaskarJoshi-01/Distributed-Systems/blob/main/Assignment2/Q1/tester.py .


# Q2
- Requirements: Find the connected components in the graph. You will be required to give the Nodes in each component in each new line. Component of size 1 is possible.
- Input: The input file will contain edges only.
- Output: The output should be a set of lines where each line contains the
nodes in a connected component.

Runner Script:
Basically we do not require any runner script , the input from should be directly redirected to mapper but something has to be done in this file before submission and according to submission's demand xD.

Mapper.py:
In this I have used DFS to get connected components by making sets if nodes are connected and if there is an intersections of sets, we just merge them and increase the value of connected component count (ccc) by one.
So what we have as nodes, we just run dfs on it. The important thing to tell is I have used defaultdict and key is node and value is its adj list.

Reducer.py:
What we obtain from mapper is node and adj list baiscally the nodes that are connected together, now what our runner will do is to take all these nodes which were the keys and find intersection of sets , if there is intersection of the sets we will add the two nodes. 
As we had recieved clarification, we have nodes values limited to 100 so in the end I have handled the case of single node by just running a loop till 100 and checking if there is some adj list associated with it and prints them in output.txt file.
Thus we will have the different connected components !


I have also made tester script for the above question. What it does it asks the number of edges we require and we can put 100C2 values in it as there are 100 nodes. Also there is constraint of 1000 edges in question idk why so simple.
So in the scirpt there is a Graph class that calculates the number of connected components(from gfg) and we just compare our ans to actual answer.







> Note : This readme was written 10 days before submission deadline so it might have some updates regarding submission commands or scripts. Anyways the Assingment can be found on https://github.com/BhaskarJoshi-01/Distributed-Systems/tree/main/Assignment2 after submission deadline or contact me on bhaskar.joshi@research.iiit.ac.in to get view permissions !



