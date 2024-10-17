#Start from bionic base Ubuntu image.
FROM ubuntu:20.04

#ADD . .

RUN apt-get update \
    && apt-get install -y python pip libatlas-base-dev \
    && apt-get clean

RUN pip install python-dotenv pybind11 scipy python-config
# run the command
