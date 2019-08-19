CXXFLAGS = -O4 -D__ZLIB_AVAILABLE__ -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS

# C++ flags for different platform
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
   #CXXFLAGS += -Werror -Wno-unused-variable
   CXXFLAGS += -Werror -Wno-unused-variable -Wno-unused-result
endif
ifeq ($(UNAME), Darwin)
   CXXFLAGS += -shared
endif

SOURCE := $(wildcard *.cpp)
OBJ = $(SOURCE:.cpp=.o)
LIB = libcore.a

all: $(LIB)

$(LIB): $(OBJ)
	ar -cr $(LIB) $(OBJ)

%.o: %.cpp
	g++ -c $(CXXFLAGS) $<
clean:
	rm -f *.o $(LIB)
