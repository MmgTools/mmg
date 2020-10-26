CC = g++  # C compiler
CFLAGS = -fPIC -g # C flags
LDFLAGS = -shared  # linking flags
RM = rm -f   # rm command
TARGET = pymmg$(shell python3-config --extension-suffix)  # target lib
SRC_DIR = ./
INC = $(shell python3 -m pybind11 --includes) -I../build/include/ -I../build/src/common
LIBS = -L../build/lib/ -lmmg3d
SRCS := $(wildcard $(SRC_DIR)/*.cpp)  # source files
OBJS = $(SRCS:.cpp=.o)  # objectives

.PHONY: all
all: ${TARGET}

# compile the sources
# must be done this way for one at a time compilation
# cannot use, $OBJS: $SRCS
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

# then link together as library remember to include ext libraries
$(TARGET): $(OBJS)
	$(CC) $(INC) ${LDFLAGS} -o $@ $^ $(LIBS)

.PHONY: clean
clean:	
	rm -f $(OBJS) $(TARGET)
