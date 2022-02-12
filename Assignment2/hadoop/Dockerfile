FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y default-jdk openssh-server openssh-client sudo vim && apt clean

RUN addgroup hadoop
RUN adduser --ingroup hadoop hdoop
RUN adduser hdoop sudo
RUN echo "hdoop ALL=(ALL:ALL) NOPASSWD: ALL" | tee /etc/sudoers.d/hdoop
USER hdoop
WORKDIR /home/hdoop

RUN ssh-keygen -t rsa -P "" -f ~/.ssh/id_rsa
RUN cat .ssh/id_rsa.pub >> .ssh/authorized_keys

RUN wget https://dlcdn.apache.org/hadoop/common/hadoop-3.3.1/hadoop-3.3.1.tar.gz \
    && tar -xf hadoop-3.3.1.tar.gz && mv hadoop-3.3.1 hadoop && rm -f hadoop-3.3.1.tar.gz

COPY conf /home/hdoop/conf
RUN (cd /home/hdoop/conf && bash configure.sh)

ENTRYPOINT [ "bash", "/home/hdoop/conf/init.sh"]
