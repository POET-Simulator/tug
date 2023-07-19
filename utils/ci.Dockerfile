FROM alpine
    MAINTAINER Max Luebke <mluebke@uni-potsdam.de>

RUN apk add --no-cache build-base openmp cmake git eigen-dev
