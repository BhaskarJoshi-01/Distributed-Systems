test: input.txt
	python3 runner.py > input
	hdfs dfs -copyFromLocal input /
	hadoop jar ../jar_files/hadoop-streaming-3.3.1.jar -input /input -output /output -mapper "python3 /mnt/Q1/mapper.py" -reducer "python3 /mnt/Q1/reducer.py"
	hdfs dfs -getmerge  /output/* /mnt/Q1/output.txt
	# diff -s /mnt/Q1/output.txt /mnt/Q1/output.txt.expected

clean:
	hdfs dfs -rm -r /input
	hdfs dfs -rm -r /output
	rm output.txt
	rm .output.txt.crc
