

CXX = mpic++
CXX_FLAGS = -std=c++11 -v -larmadillo -framework Accelerate
CXX_OPTS = -O3
BUILD_DIR = exe
EXE = qGeod
SRC = src/front/qGeod.cpp


all:
	mkdir $(BUILD_DIR)
	$(CXX) $(CXX_FLAGS) $(CXX_OPTS) $(SRC) -o $(BUILD_DIR)/$(EXE)

.PHONY: clean

clean:
	rm -rf $(BUILD_DIR)
