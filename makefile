

CXX = mpic++
CXX_FLAGS = -std=c++11 -v -larmadillo -framework Accelerate -fopenmp -D_GLIBCXX_PARALLEL
CXX_OPTS = -O3
BUILD_DIR = exe
EXE = qGeod
SRC = src/front/qGeod.cpp


all:
	mkdir $(BUILD_DIR)
	$(CXX) $(CXX_FLAGS) $(CXX_OPTS) $(SRC) -o $(BUILD_DIR)/$(EXE)
	cp input/*.mat exe/	

.PHONY: clean

clean:
	rm -rf $(BUILD_DIR)
