# (https://github.com/cni/gephysio)
#
# For parsing and filtering the physiological data (PPG 
# and respiratory data) recorded during an EPI time series.
# First build the Matlab Runtime standalone application,
# then run in command-line:
#
# docker build --no-cache -t gephysio .
#
#
# Start with the Matlab r2018b runtime container
FROM flywheel/matlab-mcr:v95
# Install Python
RUN apt-get -qq update && apt-get -qq -y install \
    python python-pip 
# for MATLAB Signal Processing Toolbox (findpeaks)
RUN apt-get install -y libasound2 libatk1.0-0
RUN pip install flywheel-sdk~=14.6

# Copy shell script and MCR binary to /bin
COPY gephysio /bin
COPY run_gephysio.sh /bin

# change permissions 
RUN chmod +x /bin/*

# Make directory for flywheel spec
ENV FLYWHEEL /flywheel/v0
RUN mkdir -p ${FLYWHEEL}

# Copy and configure run script and metadata code
COPY run.py ${FLYWHEEL}/run
RUN chmod +x ${FLYWHEEL}/run
COPY manifest.json ${FLYWHEEL}/manifest.json

# Configure entrypoint
ENTRYPOINT ["/flywheel/v0/run"]
