COMPILER = g++
CFLAGS = -O3 -mcmodel=large -std=c++11 -I./minisat/ -pthread

PRE_CFLAGS = ${CFLAGS} -c
TARGET = Generator

SRC_DIR = src

BIN_DIR = bin

EXPANDOR = Expandor
EXPANDOR_TARGET = ${SRC_DIR}/${EXPANDOR}.o
EXPANDOR_CPP_FILE = ${SRC_DIR}/${EXPANDOR}.cpp
EXPANDOR_H_FILE = ${SRC_DIR}/${EXPANDOR}.h
EXPANDOR_SOURCE_FILES = ${EXPANDOR_H_FILE} ${EXPANDOR_CPP_FILE}

OPTIMIZER = Optimizer
OPTIMIZER_TARGET = ${SRC_DIR}/${OPTIMIZER}.o
OPTIMIZER_CPP_FILE = ${SRC_DIR}/${OPTIMIZER}.cpp
OPTIMIZER_H_FILE = ${SRC_DIR}/${OPTIMIZER}.h
OPTIMIZER_SOURCE_FILES = ${OPTIMIZER_H_FILE} ${OPTIMIZER_CPP_FILE}

GLOBALOPTIMIZER = GlobalOptimizer
GLOBALOPTIMIZER_TARGET = ${SRC_DIR}/${GLOBALOPTIMIZER}.o
GLOBALOPTIMIZER_CPP_FILE = ${SRC_DIR}/${GLOBALOPTIMIZER}.cpp
GLOBALOPTIMIZER_H_FILE = ${SRC_DIR}/${GLOBALOPTIMIZER}.h
GLOBALOPTIMIZER_SOURCE_FILES = ${GLOBALOPTIMIZER_H_FILE} ${GLOBALOPTIMIZER_CPP_FILE}

EXT_MINISAT_TARGET = minisat_ext/Solver.o minisat_ext/Ext.o minisat_ext/BlackBoxSolver.o
EXT_MINISAT_SOURCE_FILES = minisat_ext/Solver.cc minisat_ext/Solver.h minisat_ext/Ext.h minisat_ext/Ext.cc minisat_ext/BlackBoxSolver.h minisat_ext/BlackBoxSolver.cc

TARGET_FILES = 	${EXPANDOR_TARGET} \
				${OPTIMIZER_TARGET} \
				${GLOBALOPTIMIZER_TARGET} \
				${EXT_MINISAT_TARGET} \

UPDATE = update
CLEAN = clean
CLEANUP = cleanup

MAIN_SOURCE_FILE = ${SRC_DIR}/main.cpp

all: ${TARGET_FILES} ${TARGET} ${UPDATE} ${CLEAN}

${EXPANDOR_TARGET}: ${EXPANDOR_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${EXPANDOR_CPP_FILE} minisat/core/Solver.o -o ${EXPANDOR_TARGET}

${OPTIMIZER_TARGET}: ${OPTIMIZER_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${OPTIMIZER_CPP_FILE} minisat/core/Solver.o -o ${OPTIMIZER_TARGET}

${GLOBALOPTIMIZER_TARGET}: ${GLOBALOPTIMIZER_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${GLOBALOPTIMIZER_CPP_FILE} minisat/core/Solver.o -o ${GLOBALOPTIMIZER_TARGET}

${EXT_MINISAT_TARGET}: ${EXT_MINISAT_SOURCE_FILES}
	cd minisat_ext; make -f Makefile extminisat; cd ..

${TARGET}: ${MAIN_SOURCE_FILE} ${TARGET_FILES}
	${COMPILER} ${CFLAGS} ${MAIN_SOURCE_FILE} ${TARGET_FILES} -o ${TARGET}

${UPDATE}:
	chmod +x ${BIN_DIR}/*

${UPDATE}:
	chmod +x ${BIN_DIR}/*

${CLEAN}:
	rm -f *~
	rm -f ${SRC_DIR}/*.o
	rm -f ${SRC_DIR}/*~
	rm -f minisat/utils/*.or minisat/utils/*.o minisat/core/*.or minisat/core/*.o
	rm -f minisat/core/depend.mk
	rm -f ${EXT_MINISAT_TARGET}
	rm -f minisat_ext/depend.mk

${CLEANUP}: ${CLEAN}
	rm -f ${TARGET}