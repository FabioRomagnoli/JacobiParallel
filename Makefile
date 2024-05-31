MAKEFILEH_DIR=../
include $(MAKEFILEH_DIR)/Makefile.inc
-include Makefile.inc


CXX = mpic++
CXXFLAGS += -fopenmp
INCLUDES += -I$(mkEigenInc) -I$(MUPARSER_INCLUDE)  #-I$(UTIL_INCLUDE)


#REQUIRMENTS
LIBD+=$(LDLIBS)

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
	$(CXX)  $(CXXFLAGS)  $(CPPFLAGS) $(INCLUDES) -o $@ $^ $(LIBD)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) $(EXEC)
	$(RM) *.o

distclean: clean
	$(RM) $(EXEC)