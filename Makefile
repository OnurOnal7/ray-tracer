CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O3 -march=native -pthread

TARGET = main

SRC = main.cpp aabb.cpp bvh.cpp loader.cpp
OBJS = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run
