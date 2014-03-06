CC=g++
SOURCES=kernelTest.cpp SMOClassifierTest.cpp ModifiedSMOClassifierTest.cpp AllTests.cpp AsmFindHighLowTest.cpp
OBJECTS=$(patsubst %.cpp,%.o,$(SOURCES))
LDFLAGS=-ldl -pthread
CFLAGS=-g -m32 -Wall -Wextra -std=c++11
INCLUDE=-I../contrib -I../src
OUTPUT=tests
GTEST_DIR=../contrib
GTEST_SOURCE=gtest/gtest-all.cc
GTEST_LINK=gtest-all.o
all: $(GTEST_LINK) $(OUTPUT) 

$(OUTPUT) : $(OBJECTS)
	$(CC) -g -m32 -Wall $(INCLUDE)  -o $(OUTPUT) $(OBJECTS) $(GTEST_LINK) $(LDFLAGS)
	
$(GTEST_LINK) :
	$(CC) $(CFLAGS) $(INCLUDE) -c $(GTEST_DIR)/$(GTEST_SOURCE) -o $@ -lpthread

.cpp.o:
	$(CC)  $(CFLAGS) $(INCLUDE) -c $< -o $@

clean:
	rm -rf $(OBJECTS);
