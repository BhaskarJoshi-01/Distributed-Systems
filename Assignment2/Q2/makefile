test: input
	hdfs dfs -copyFromLocal input /
	hadoop jar ../jar_files/hadoop-streaming-3.3.1.jar -input /input -output /output -mapper "python3 /mnt/Q2/mapper.py" -reducer "python3 /mnt/Q2/reducer.py"

clean:
	hdfs dfs -rm -r /input
	hdfs dfs -rm -r /output
