#FROM ubuntu:20.04
FROM  public.ecr.aws/lambda/python:3.9

#ARG DEBIAN_FRONTEND=noninteractive
#RUN apt-get update && \
#    apt-get install -y --fix-missing dirmngr gnupg apt-transport-https ca-certificates software-properties-common
#RUN apt-get install -y --fix-missing git make cmake python3 python3-pip maven tabix gcc apt-transport-https curl openssl libssl-dev zip

RUN yum -y update && \
    yum -y groupinstall "Development Tools" && \
    yum -y install \
    epel-release \
    chkconfig \
    alternatives \
    dirmngr \
    gnupg \
    ca-certificates \
    make \
    cmake3 \
    maven \
    gcc \
    gcc-c++ \
    curl \
    libcurl \
    libcurl-devel \
    wget \
    openssl \
    dnf-plugins-core \
    yum-utils \
    zip \
    gzip \
    bzip2 \
    bzip2-devel \
    tar \
    xz \
    xz-devel \
    zlib-devel \
    zlib && \
    yum clean all

#Tabix & Bgzip
ENV HTSLIB_VERSION=1.12
ADD https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 /tmp/
RUN tar -xvjf /tmp/htslib-${HTSLIB_VERSION}.tar.bz2 -C /tmp && \
    cd /tmp/htslib-${HTSLIB_VERSION} && \
    make && \
    make install && \
    cd .. && \
    rm -rf /tmp/htslib-${HTSLIB_VERSION} /tmp/htslib-${HTSLIB_VERSION}.tar.bz2


RUN ln -sf /usr/bin/cmake3 /usr/bin/cmake
RUN pip install --upgrade cget pandas certifi pysam pybase64 jsons boto3 regex pip-system-certs statsmodels scipy
#RUN pip3 install --upgrade pip
#RUN pip3 install --upgrade certifi
#RUN python3 --version && python3 -m certifi

#Add python scripts
COPY lambda.py ${LAMBDA_TASK_ROOT}
COPY file_io.py ${LAMBDA_TASK_ROOT}
COPY impute.py ${LAMBDA_TASK_ROOT}
COPY prs.py ${LAMBDA_TASK_ROOT}
COPY logging_config.py ${LAMBDA_TASK_ROOT}

#Build Minimac4 from GitHub
RUN cget install --insecure --prefix /usr/local/minimac4 statgen/Minimac4@v4.1.6
RUN chmod +x /usr/local/minimac4/bin/minimac4

#Build Eagle2
ADD https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz /usr/local/
RUN cd /usr/local/ && \
    tar xzvf Eagle_v2.4.1.tar.gz && \
    chmod +x /usr/local/Eagle_v2.4.1/eagle

    #Add minimac4 and eagle to path
ENV PATH="${PATH}:/usr/local/minimac4/bin:/usr/local/Eagle_v2.4.1:/usr/local/bin"

# Call the AWS Lambda handler function
CMD ["lambda.handler"]
