FROM python:3.12-slim

RUN apt update && apt upgrade -y && apt install -y procps

WORKDIR /hairpin2

COPY . /hairpin2

RUN pip install --root-user-action ignore --no-warn-script-location /hairpin2

RUN LOC=$(which hairpin2) \
    && if [ -z "$LOC" ]; then \
    echo "hairpin install failed" && exit 1; \
    else echo "hairpin install successful"; fi

# Set up the default command for the container
ENTRYPOINT ["hairpin2"]

