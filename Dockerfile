FROM python:3.12-slim

RUN apt update && apt upgrade -y && apt install -y procps

WORKDIR /hairpin2

COPY . /hairpin2

RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu

RUN pip install --no-warn-script-location /hairpin2

ENV PATH=$PATH:/home/ubuntu/.local/bin

# Define a test script to check the installation of hairpin
RUN LOC=$(which hairpin2) \
    && if [ -z "$LOC" ]; then \
    echo "hairpin install failed" && exit 1; \
    else echo "hairpin install successful"; fi

# Set up the default command for the container
ENTRYPOINT ["hairpin2"]

