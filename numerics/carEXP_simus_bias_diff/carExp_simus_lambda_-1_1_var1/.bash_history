ls
tar -xzf R402.tar.gz
export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R
R --version
mkdir packages
R_LIBS=$PWD/packages
R
ls
R
tar -czf packages.tar.gz packages/
ls
exit
