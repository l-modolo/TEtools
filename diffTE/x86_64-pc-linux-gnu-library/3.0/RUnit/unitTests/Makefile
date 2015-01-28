TOP = ../../..
TEST_DIR = `pwd`
INST_PATH = `dirname ${TEST_DIR}`
PKG_PATH = `dirname ${INST_PATH}`

## Note that 'R' was not right here: there may be no 'R' on the path.
all: install test

install:
	cd ${TOP}; \
	$(R_HOME)/bin/R CMD INSTALL --clean ${PKG_PATH} && \
	cd ${TEST_DIR}

test:
	$(R_HOME)/bin/R --slave < runalltests.R
