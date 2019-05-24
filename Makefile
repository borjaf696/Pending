CXXFLAGS = -Wall -std=c++14 -g -O3
LBOOST = -lz -lboost_system -lboost_filesystem -lboost_regex -lboost_program_options
LDSDSL = -lsdsl -ldivsufsort -ldivsufsort64
LINCLUDE = -I/home/bfreire/esprela/sdsl/include/

BIN_release = bin/esprela

#Print_Execution
p_obj := ${patsubst %.cpp,%.o,${wildcard utils/*.cpp}}
utils/%.o: utils/%.cpp utils/*h
	${CXX} -c ${CXXFLAGS} ${LBOOST} $< -o $@

#Main_Execution
m_main_obj := ${patsubst %.cpp,%.o,${wildcard src/*.cpp}}
src/%.o: src/%.cpp utils/*h
	${CXX} -c ${CXXFLAGS} ${LINCLUDE} ${LDSDSL} ${LBOOST} $< -o $@
m_code := src/main.o

#Running
all: ${p_obj} ${m_main_obj}
	${CXX} ${LINCLUDE} ${LBOOST} ${LDSDSL} ${p_obj} ${m_code} -o ${BIN_release}

#Clean
clean:
	-rm ${p_obj}
	-rm ${m_main_obj}
	-rm ${BIN_release}
