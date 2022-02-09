cat .bashrc >> ~/.bashrc
source .bashrc
echo 'export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64' >> $HADOOP_HOME/etc/hadoop/hadoop-env.sh
sudo mv hadoop-env.sh \
    core-site.xml \
    hdfs-site.xml \
    mapred-site.xml \
    yarn-site.xml $HADOOP_HOME/etc/hadoop/

hdfs namenode -format
