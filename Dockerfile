FROM public.ecr.aws/lambda/python:3.9

# Set up basic dev tools and system dependencies early (rarely changes)
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
    curl-devel \
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
    bcftools \
    htslib \
    tabix \
    tar \
    which \
    xz \
    xz-devel \
    zlib-devel \
    zlib && \
    yum clean all && \
    ln -sf /usr/bin/cmake3 /usr/bin/cmake

# ---------------------------
# Install Python packages
# ---------------------------
COPY requirements.txt .
RUN pip install --upgrade pip && \
    pip install -r requirements.txt

# ---------------------------
# Install HTSLIB from source
# ---------------------------
RUN curl -L -O https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 \
 && tar -xvjf htslib-1.19.1.tar.bz2 \
 && cd htslib-1.19.1 && ./configure && make && make install

# ---------------------------
# Build Minimac4 from GitHub
# ---------------------------
RUN cget install --insecure --prefix /usr/local/minimac4 statgen/Minimac4@v4.1.6 && \
    chmod +x /usr/local/minimac4/bin/minimac4

# ---------------------------
# Download and extract Eagle2
# ---------------------------
ADD https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz /usr/local/
RUN cd /usr/local/ && \
    tar xzvf Eagle_v2.4.1.tar.gz && \
    chmod +x /usr/local/Eagle_v2.4.1/eagle

# ---------------------------
# Build BCFtools from source
# ---------------------------
RUN cd /usr/local && \
    wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 && \
    tar -xvjf bcftools-1.22.tar.bz2 && \
    cd bcftools-1.22/ && \
    make && make install && \
    chmod +x /usr/local/bcftools-1.22/bcftools

# ---------------------------
# Add Eagle and Minimac4 to PATH
# ---------------------------
ENV PATH="${PATH}:/usr/local/minimac4/bin:/usr/local/Eagle_v2.4.1:/usr/local/bin"

# ---------------------------
# Copy your function code last (this changes most often)
# ---------------------------
COPY lambda.py ${LAMBDA_TASK_ROOT}
COPY file_io.py ${LAMBDA_TASK_ROOT}
COPY impute.py ${LAMBDA_TASK_ROOT}
COPY prs.py ${LAMBDA_TASK_ROOT}
COPY logging_config.py ${LAMBDA_TASK_ROOT}

# Default CMD to call your Lambda handler
CMD ["lambda.handler"]