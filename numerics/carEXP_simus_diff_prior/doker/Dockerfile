FROM rocker/r-ver:4.0
RUN install2.r --error \
    Rcpp \
    RcppArmadillo \
    MASS \
    Matrix \
    CholWishart \
    FNN \
    coda \ 
    igraph \
    ggraph \
    ggplot2 \
    RcppProgress

COPY CARlasso_0.1.2.tar.gz /tmp
RUN apt update
RUN apt install -y libxml2
RUN apt install -y libmkl-rt
RUN apt install -y libglpk-dev
RUN install2.r --error /tmp/CARlasso_0.1.2.tar.gz


