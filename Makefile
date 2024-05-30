CXX      ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -Wno-deprecated-enum-enum-conversion
INCLUDES ?= -I$(mkEigenInc)  #-I$(include folders) or mk modules

# MAKEFILEH_DIR=.
# include $(MAKEFILEH_DIR)/Makefile.inc
# -include Makefile.inc

# MP_PRFX=$(PACS_ROOT)
# MP_LIB_DIR=$(MP_PRFX)/lib
# MP_INCLD_DIR=$(MP_PRFX)/include

#REQUIRMENTS
# MP_REQ_LIBS = #-L$(MP_LIB_DIR) -l any lib that is needed 

.PHONY : all $(EXEC) $(OBJS) clean distclean 

# get all files *.cpp
SRCS=$(wildcard *.cpp)
# get the corresponding object file
OBJS = $(SRCS:.cpp=.o)
# get all headers in the working directory
HEADERS=$(wildcard *.hpp)

exe_sources=$(filter main%.cpp,$(SRCS))
EXEC=$(exe_sources:.cpp=)

all: $(EXEC) 

$(EXEC): $(OBJS)
	$(CXX)  $(CXXFLAGS)  $(CPPFLAGS) $(inputs) -o $@ $^ $(MP_REQ_LIBS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) $(EXEC)
	$(RM) *.o

distclean: clean
	$(RM) $(EXEC)