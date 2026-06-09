CXX ?= c++
CXXFLAGS ?= -O2
CXXFLAGS += -std=c++17 -I. -I./json/include -Wall -Wextra -Wpedantic

TARGET := ccm_truth_rates
SOURCES := ccm_truth_rates.cc FormFactor.cc NuFlux.cc xscns.cc
OBJECTS := $(SOURCES:.cc=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	$(RM) $(TARGET) $(OBJECTS)
