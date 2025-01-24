COMPILER = g++
CFLAGS = -DNDEBUG -O3 -mcmodel=large -std=c++11 -pthread -I./minisat/

PRE_CFLAGS = ${CFLAGS} -c
TARGET = Optimizer

SRC_DIR = src

SAT = SAT
SAT_TARGET = ${SRC_DIR}/${SAT}.o
SAT_CPP_FILE = ${SRC_DIR}/${SAT}.cc
SAT_H_FILE = ${SRC_DIR}/${SAT}.h
SAT_SOURCE_FILES = ${SAT_H_FILE} ${SAT_CPP_FILE}

# ConstraintFile
ConstraintFile = ConstraintFile
ConstraintFile_TARGET = ${SRC_DIR}/${ConstraintFile}.o
ConstraintFile_CPP_FILE = ${SRC_DIR}/${ConstraintFile}.C
ConstraintFile_H_FILE = ${SRC_DIR}/${ConstraintFile}.H
ConstraintFile_SOURCE_FILES = ${ConstraintFile_H_FILE} ${ConstraintFile_CPP_FILE}

# Options.o
Options = Options
Options_TARGET = ${SRC_DIR}/${Options}.o
Options_CPP_FILE = ${SRC_DIR}/${Options}.cc
Options_H_FILE = ${SRC_DIR}/${Options}.h
Options_SOURCE_FILES = ${Options_H_FILE} ${Options_CPP_FILE}

# SpecificationFile.o
SpecificationFile = SpecificationFile
SpecificationFile_TARGET = ${SRC_DIR}/${SpecificationFile}.o
SpecificationFile_CPP_FILE = ${SRC_DIR}/${SpecificationFile}.cc
SpecificationFile_H_FILE = ${SRC_DIR}/${SpecificationFile}.h
SpecificationFile_SOURCE_FILES = ${SpecificationFile_H_FILE} ${SpecificationFile_CPP_FILE}

# PascalTriangle.o
PascalTriangle = PascalTriangle
PascalTriangle_TARGET = ${SRC_DIR}/${PascalTriangle}.o
PascalTriangle_CPP_FILE = ${SRC_DIR}/${PascalTriangle}.cc
PascalTriangle_H_FILE = ${SRC_DIR}/${PascalTriangle}.h
PascalTriangle_SOURCE_FILES = ${PascalTriangle_H_FILE} ${PascalTriangle_CPP_FILE}

# Combinadic.o
Combinadic = Combinadic
Combinadic_TARGET = ${SRC_DIR}/${Combinadic}.o
Combinadic_CPP_FILE = ${SRC_DIR}/${Combinadic}.cc
Combinadic_H_FILE = ${SRC_DIR}/${Combinadic}.h
Combinadic_SOURCE_FILES = ${Combinadic_H_FILE} ${Combinadic_CPP_FILE}

# LocalSearch.o
LocalSearch = LocalSearch
LocalSearch_TARGET = ${SRC_DIR}/${LocalSearch}.o
LocalSearch_CPP_FILE = ${SRC_DIR}/${LocalSearch}.cc
LocalSearch_H_FILE = ${SRC_DIR}/${LocalSearch}.h
LocalSearch_SOURCE_FILES = ${LocalSearch_H_FILE} ${LocalSearch_CPP_FILE}

# Coverage.o
Coverage = Coverage
Coverage_TARGET = ${SRC_DIR}/${Coverage}.o
Coverage_CPP_FILE = ${SRC_DIR}/${Coverage}.cc
Coverage_H_FILE = ${SRC_DIR}/${Coverage}.h
Coverage_SOURCE_FILES = ${Coverage_H_FILE} ${Coverage_CPP_FILE}

# CoveringArray.o
CoveringArray = CoveringArray
CoveringArray_TARGET = ${SRC_DIR}/${CoveringArray}.o
CoveringArray_CPP_FILE = ${SRC_DIR}/${CoveringArray}.cc
CoveringArray_H_FILE = ${SRC_DIR}/${CoveringArray}.h
CoveringArray_SOURCE_FILES = ${CoveringArray_H_FILE} ${CoveringArray_CPP_FILE}

# TupleSet.o
TupleSet = TupleSet
TupleSet_TARGET = ${SRC_DIR}/${TupleSet}.o
TupleSet_CPP_FILE = ${SRC_DIR}/${TupleSet}.cc
TupleSet_H_FILE = ${SRC_DIR}/${TupleSet}.h
TupleSet_SOURCE_FILES = ${TupleSet_H_FILE} ${TupleSet_CPP_FILE}

# OptionTupleSet.o
OptionTupleSet = OptionTupleSet
OptionTupleSet_TARGET = ${SRC_DIR}/${OptionTupleSet}.o
OptionTupleSet_CPP_FILE = ${SRC_DIR}/${OptionTupleSet}.cc
OptionTupleSet_H_FILE = ${SRC_DIR}/${OptionTupleSet}.h
OptionTupleSet_SOURCE_FILES = ${OptionTupleSet_H_FILE} ${OptionTupleSet_CPP_FILE}

# LineVarTupleSet.o
LineVarTupleSet = LineVarTupleSet
LineVarTupleSet_TARGET = ${SRC_DIR}/${LineVarTupleSet}.o
LineVarTupleSet_CPP_FILE = ${SRC_DIR}/${LineVarTupleSet}.cc
LineVarTupleSet_H_FILE = ${SRC_DIR}/${LineVarTupleSet}.h
LineVarTupleSet_SOURCE_FILES = ${LineVarTupleSet_H_FILE} ${LineVarTupleSet_CPP_FILE}

# mersenne.o
mersenne = mersenne
mersenne_TARGET = ${SRC_DIR}/${mersenne}.o
mersenne_CPP_FILE = ${SRC_DIR}/${mersenne}.cc
mersenne_H_FILE = ${SRC_DIR}/${mersenne}.h
mersenne_SOURCE_FILES = ${mersenne_H_FILE} ${mersenne_CPP_FILE}

EXT_MINISAT_TARGET = minisat_ext/Solver.o minisat_ext/Ext.o minisat_ext/BlackBoxSolver.o
EXT_MINISAT_SOURCE_FILES = minisat_ext/Solver.cc minisat_ext/Solver.h minisat_ext/Ext.h minisat_ext/Ext.cc minisat_ext/BlackBoxSolver.h minisat_ext/BlackBoxSolver.cc

TARGET_FILES = 	${SAT_TARGET} \
				${ConstraintFile_TARGET} \
				${Options_TARGET} \
				${SpecificationFile_TARGET} \
				${PascalTriangle_TARGET} \
				${Combinadic_TARGET} \
				${LocalSearch_TARGET} \
				${Coverage_TARGET} \
				${CoveringArray_TARGET} \
				${TupleSet_TARGET} \
				${OptionTupleSet_TARGET} \
				${LineVarTupleSet_TARGET} \
				${mersenne_TARGET} \
				${EXT_MINISAT_TARGET} \

CLEAN = clean
CLEANUP = cleanup

MAIN_SOURCE_FILE = ${SRC_DIR}/main.cc

all: ${TARGET_FILES} ${TARGET} ${CLEAN}

${EXT_MINISAT_TARGET}: ${EXT_MINISAT_SOURCE_FILES}
	cd minisat_ext; make -f Makefile extminisat; cd ..

${SAT_TARGET}: ${SAT_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${SAT_CPP_FILE} minisat/core/Solver.o -o ${SAT_TARGET}

${ConstraintFile_TARGET}: ${ConstraintFile_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${ConstraintFile_CPP_FILE} minisat/core/Solver.o -o ${ConstraintFile_TARGET}

${Options_TARGET}: ${Options_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${Options_CPP_FILE} -o ${Options_TARGET}

${SpecificationFile_TARGET}: ${SpecificationFile_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${SpecificationFile_CPP_FILE} -o ${SpecificationFile_TARGET}

${PascalTriangle_TARGET}: ${PascalTriangle_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${PascalTriangle_CPP_FILE} -o ${PascalTriangle_TARGET}

${Combinadic_TARGET}: ${Combinadic_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${Combinadic_CPP_FILE} -o ${Combinadic_TARGET}

${LocalSearch_TARGET}: ${LocalSearch_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${LocalSearch_CPP_FILE} -o ${LocalSearch_TARGET}

${Coverage_TARGET}: ${Coverage_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${Coverage_CPP_FILE} minisat/core/Solver.o -o ${Coverage_TARGET}

${CoveringArray_TARGET}: ${CoveringArray_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${CoveringArray_CPP_FILE} minisat/core/Solver.o -o ${CoveringArray_TARGET}

${TupleSet_TARGET}: ${TupleSet_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${TupleSet_CPP_FILE} -o ${TupleSet_TARGET}

${OptionTupleSet_TARGET}: ${OptionTupleSet_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${OptionTupleSet_CPP_FILE} -o ${OptionTupleSet_TARGET}

${LineVarTupleSet_TARGET}: ${LineVarTupleSet_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${LineVarTupleSet_CPP_FILE} -o ${LineVarTupleSet_TARGET}

${mersenne_TARGET}: ${mersenne_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${mersenne_CPP_FILE} -o ${mersenne_TARGET}

${TARGET}: ${MAIN_SOURCE_FILE} ${TARGET_FILES}
	${COMPILER} ${CFLAGS} ${MAIN_SOURCE_FILE} ${TARGET_FILES} -o ${TARGET}

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

