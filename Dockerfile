FROM rocker/r-ver:4.3.3

RUN apt-get update && \
    apt-get install -y libgoogle-perftools4 && \
    rm -rf /var/lib/apt/lists/*
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libtcmalloc.so.4

RUN install2.r --error plumber caret dplyr pROC ranger jsonlite

WORKDIR /app

COPY cancer_api.R /app/
COPY models/ /app/models/

EXPOSE 8000
CMD ["R", "-e", "pr <- plumber::plumb('cancer_api.R'); pr$run(host='0.0.0.0', port=8000)"]