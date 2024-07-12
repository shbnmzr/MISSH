# Define base directories and compiler
DEBUG_DIR := ./debug
BUILD_DIR := ./build
SOURCE_DIR := ./src
GPP_PATH := /opt/homebrew/bin/g++-14
CC := $(GPP_PATH)
CXX := $(GPP_PATH)

# If make debug was called, define directories accordingly and add -g flag
ifeq (debug,$(filter debug,$(MAKECMDGOALS)))
    OBJ_DIR := $(DEBUG_DIR)/obj
    DEPS_DIR := $(DEBUG_DIR)/deps

    CFLAGS = -g -std=c++0x -O3 -Wall -c -fmessage-length=0 -MMD -MP -fopenmp

# If make clean all, define directories accordingly
else ifeq (all,$(filter all,$(MAKECMDGOALS)))
    OBJ_DIR := $(BUILD_DIR)/obj
    DEPS_DIR := $(BUILD_DIR)/deps

    CFLAGS = -std=c++0x -O3 -Wall -c -fmessage-length=0 -MMD -MP -fopenmp
else
    OBJ_DIR := $(BUILD_DIR)/obj
    DEPS_DIR := $(BUILD_DIR)/deps

    CFLAGS = -std=c++0x -O3 -Wall -c -fmessage-length=0 -MMD -MP -fopenmp
endif

# Ensure OBJ_DIR and DEPS_DIR are always defined
OBJ_DIR ?= $(BUILD_DIR)/obj
DEPS_DIR ?= $(BUILD_DIR)/deps

# CPLEX library paths
CPLEX_DIR := /Applications/CPLEX_Studio2211
CPLEX_INCLUDE := $(CPLEX_DIR)/cplex/include
CONCERT_INCLUDE := $(CPLEX_DIR)/concert/include
CPLEX_LIB := $(CPLEX_DIR)/cplex/lib/arm64_osx/static_pic
CONCERT_LIB := $(CPLEX_DIR)/concert/lib/arm64_osx/static_pic

INCLUDES := -I$(CPLEX_INCLUDE) -I$(CONCERT_INCLUDE)
LIBS := -L$(CPLEX_LIB) -L$(CONCERT_LIB) -lcplex -lilocplex -lconcert -lpthread

# If on unix, search for all .cpp files from sources directory
CPP_SRCS := $(shell find $(SOURCE_DIR) -name "*.cpp")

# Create objs and deps lists, and a list with all subdirectories to create inside OBJ_DIR
OBJS := $(CPP_SRCS:$(SOURCE_DIR)/%.cpp=$(OBJ_DIR)/%.o)
CPP_DEPS := $(CPP_SRCS:$(SOURCE_DIR)/%.cpp=$(OBJ_DIR)/%.d)
TREE_UNIX := $(sort $(patsubst %/,%,$(dir $(OBJS))))

# Debugging help
$(info DEBUG_DIR=$(DEBUG_DIR))
$(info BUILD_DIR=$(BUILD_DIR))
$(info SOURCE_DIR=$(SOURCE_DIR))
$(info OBJ_DIR=$(OBJ_DIR))
$(info DEPS_DIR=$(DEPS_DIR))
$(info CPP_SRCS=$(CPP_SRCS))
$(info OBJS=$(OBJS))
$(info CPP_DEPS=$(CPP_DEPS))
$(info TREE_UNIX=$(TREE_UNIX))

# All target
all: $(BUILD_DIR)/ISSH

$(BUILD_DIR)/ISSH: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CXX) -o "$@" $(OBJS) $(USER_OBJS) $(LIBS) -fopenmp
	@echo 'Finished building target: $@'
	@echo ' '

# Debug target
debug: $(DEBUG_DIR)/ISSH

$(DEBUG_DIR)/ISSH: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CXX) -o "$@" $(OBJS) $(USER_OBJS) $(LIBS) -fopenmp
	@echo 'Finished building target: $@'
	@echo ' '

# Compile all dependencies for ISSH
.SECONDEXPANSION:
$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $$(@D)
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) $(INCLUDES) -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Create subdirectories tree, depending on which system we are
$(TREE_UNIX): %:
	mkdir -p $@
	mkdir -p $(@:$(OBJ_DIR)%=$(DEPS_DIR)%)

remake: clean-build all

remake-debug: clean-debug debug

clean: clean-build clean-debug

clean-build:
	$(RM) -r $(BUILD_DIR)
	-@echo ' '

clean-debug:
	$(RM) -r $(DEBUG_DIR)
	-@echo ' '
