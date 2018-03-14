# Makefile for the Cahn-Hillard differential equation solver.

SRC_DIR=src
HEADERS=$(wildcard $(SRC_DIR)/*.hpp)
SRC_FILES=$(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES=$(patsubst $(SRC_DIR)/%.cpp, %.o, $(SRC_FILES))


CXX=g++
CPPSTD=-std=c++14 
DEBUG=-g
OPT=-O2
LFLAGS= -lboost_program_options -lboost_system -lboost_filesystem
INC=-I$(SRC_DIR) -I$(TEST_DIR) -I$(HOME)/include

EXE_FILE=cahnHilliard



$(EXE_FILE): $(OBJ_FILES) 
	$(CXX) $(CPPSTD) $(OPT) -o $@  $^ $(LFLAGS)


## objs      : create object files
.PHONY : objs
objs : $(OBJ_FILES) $(TEST_OBJ_FILES)

%.o : $(SRC_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CPPSTD) $(OPT) -c $< -o $@ $(INC) 



## clean     : remove auto generated files
.PHONY : clean
clean :
	rm -f $(OBJ_FILES)
	rm -f $(EXE_FILE)
	rm -f *.log

## variables : Print variables
.PHONY :variables
variables:
	@echo SRC_DIR:        $(SRC_DIR)
	@echo SRC_FILES:      $(SRC_FILES)
	@echo OBJ_FILES:      $(OBJ_FILES)



## help      : Print help
.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<
