sudo /etc/init.d/ssh start
cd /home/hdoop/hadoop/sbin
./start-dfs.sh
./start-yarn.sh
jps
tail -f /dev/null
