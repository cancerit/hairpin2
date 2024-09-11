FROM python:3.12-slim

# Set the working directory inside the container
WORKDIR /hairpin2

# Copy the current working directory contents into the container
COPY . /hairpin2

# Install the hairpin package
RUN pip install --root-user-action ignore /hairpin2/

# Define a test script to check the installation of hairpin
RUN LOC=$(which hairpin2) \
    && if [ -z "$LOC" ]; then \
    echo "hairpin install failed" && exit 1; \
    else echo "hairpin install successful"; fi

# Set up the default command for the container
ENTRYPOINT ["hairpin2"]

