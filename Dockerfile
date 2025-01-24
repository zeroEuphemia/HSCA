FROM ubuntu:22.04 AS build-hsca
RUN apt-get update
RUN apt-get install -y gcc g++ python3 make
RUN apt-get install -y libz-dev
COPY . /hsca
WORKDIR /hsca
RUN sh build.sh

FROM python:3.12-slim AS release-hsca
LABEL org.opencontainers.image.source=https://github.com/zeroEuphemia/HSCA
WORKDIR /hsca
COPY benchmarks /hsca/benchmarks
COPY --from=build-hsca \
    /hsca/SamplingCA /hsca/Generator /hsca/Optimizer /hsca/run.py \
    /hsca/formatencoding.py /hsca/format_converter.sh \
    /hsca/
COPY --from=build-hsca \
    /hsca/format_converter/*.py /hsca/format_converter/FormatConverter /hsca/format_converter/ctw_parser \
    /hsca/format_converter/
COPY --from=build-hsca \
    /hsca/bin/* \
    /hsca/bin/
RUN mkdir -p tmp format_converter/tmp
ENTRYPOINT [ "python3", "run.py" ]
