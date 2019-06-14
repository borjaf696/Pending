CXXFLAGS = -Wall -std=c++14 -g -O3
LBOOST = -lz -lboost_system -lboost_filesystem -lboost_regex -lboost_program_options
LDSDSL = -lsdsl -ldivsufsort -ldivsufsort64
LINCLUDE = -I/home/bfreire/sdsl/include/

BIN_release = bin/esprela

#Print_Execution
p_obj := ${patsubst %.cpp,%.o,${wildcard utils/*.cpp}}
utils/%.o: utils/%.cpp utils/*hpp
	${CXX} -c ${CXXFLAGS} ${LBOOST} $< -o $@

#FS
fs_obj := ${patsubst %.cpp,%.o,${wildcard fs/*.cpp}}
fs/%.o: fs/%.cpp fs/*h
	${CXX} -c ${CXXFLAGS} $< -o $@

#Operations
op_obj := ${patsubst %.cpp,%.o,${wildcard operations/*.cpp}}
operations/%.o: operations/%.cpp operations/*hpp
	${CXX} -c ${CXXFLAGS} $< -o $@

#Read, seq, kmer
seq_obj := ${patsubst %.cpp,%.o,${wildcard sequence/*.cpp}}
ReadData/%.o: sequence/%.cpp sequence/*.h
	${CXX} -c ${CXXFLAGS} $< -o $@

#Main_Execution
m_main_obj := ${patsubst %.cpp,%.o,${wildcard src/*.cpp}}
src/%.o: src/%.cpp utils/*.hpp operations/*.hpp
	${CXX} -c ${CXXFLAGS} ${LINCLUDE} ${LDSDSL} ${LBOOST} $< -o $@
m_code := src/main.o

#Running
all: ${op_obj} ${p_obj} ${seq_obj} ${fs_obj} ${m_main_obj}
	${CXX} ${LINCLUDE} ${LDSDSL} ${LBOOST} ${op_obj} ${p_obj} ${seq_obj} ${fs_obj} ${m_code} -o ${BIN_release}

#Clean
clean:
	-rm ${p_obj}
	-rm ${fs_obj}
	-rm ${op_obj}
	-rm ${seq_obj}
	-rm ${m_main_obj}
	-rm ${BIN_release}
