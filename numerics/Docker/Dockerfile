FROM rocker/r-ver:4.0
RUN install2.r --error \
    Rcpp \
    RcppArmadillo \
    MASS \
    Matrix \
    CholWishart

COPY rMGIG_0.1.tar.gz /tmp

RUN install2.r --error /tmp/rMGIG_0.1.tar.gz

