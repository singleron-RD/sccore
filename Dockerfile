FROM python:3.11-slim

LABEL author="Yiqi Zhou" \
      description="single cell core" 

# Add the sccore source files to the container
ADD . /usr/src/sccore
WORKDIR /usr/src/sccore

# - Install `ps` for Nextflow
# - Install sccore through pip
# - Delete unnecessary Python files
# - Remove sccore source directory
# - Add custom group and user
RUN \
    echo "Docker build log: Run apt-get update" 1>&2 && \
    apt-get update -y -qq \
    && \
    echo "Docker build log: Install procps" 1>&2 && \
    apt-get install -y -qq procps && \
    echo "Docker build log: Clean apt cache" 1>&2 && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean -y && \
    echo "Docker build log: Upgrade pip and install sccore" 1>&2 && \
    pip install --quiet --upgrade pip && \
    pip install --verbose --no-cache-dir . && \
    echo "Docker build log: Delete python cache directories" 1>&2 && \
    find /usr/local/lib/python3.11 \( -iname '*.c' -o -iname '*.pxd' -o -iname '*.pyd' -o -iname '__pycache__' \) -printf "\"%p\" " | \
    xargs rm -rf {} && \
    echo "Docker build log: Delete /usr/src/sccore" 1>&2 && \
    rm -rf "/usr/src/sccore/" && \
    echo "Docker build log: Add sccore user and group" 1>&2 && \
    groupadd --gid 1000 sccore && \
    useradd -ms /bin/bash --create-home --gid sccore --uid 1000 sccore

# Set to be the new user
USER sccore

# Set default workdir to user home
WORKDIR /home/sccore
