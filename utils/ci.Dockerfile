FROM alpine

LABEL maintainer="Max Luebke <mluebke@uni-potsdam.de"

RUN apk add --no-cache build-base openmp cmake git eigen-dev

RUN git clone https://github.com/doctest/doctest.git /doctest \
    && cd /doctest \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make install \
    && cd / \
    && rm -rf /doctest
