LIBS = -lgsl -lgslcblas -lm
adptSTS: adptSTS.cpp latestDepSearch.cpp refPathErr.cpp HeapContainer.cpp
	${CXX} -g -std=c++11 adptSTS.cpp latestDepSearch.cpp refPathErr.cpp HeapContainer.cpp -o adptSTS ${LIBS}
