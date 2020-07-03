USER_OBJS :=
LIBS :=

# define base directories and compiler
DEBUG_DIR := ./debug
BUILD_DIR := ./build
SOURCE_DIR := ./src
CC = g++


# if make debug was called, define directories accordingly and add -g flag
ifeq (debug,$(filter debug,$(MAKECMDGOALS)))
	OBJ_DIR := $(DEBUG_DIR)/obj
	DEPS_DIR := $(DEBUG_DIR)/deps

	CFLAGS = -g -std=c++0x -O3 -Wall -c -fmessage-length=0 -MMD -MP -fopenmp

# if make clean all, define directories accordingly
else ifeq (all,$(filter all,$(MAKECMDGOALS)))
	OBJ_DIR := $(BUILD_DIR)/obj
	DEPS_DIR := $(BUILD_DIR)/deps

	CFLAGS = -std=c++0x -O3 -Wall -c -fmessage-length=0 -MMD -MP -fopenmp
endif


ifeq ($(OS),Windows_NT)
	# if on windows, search for all .cpp files from sources directory
	CPP_SRCS := $(shell FORFILES /P $(SOURCE_DIR) /S /M *.cpp /C "CMD /C ECHO @relpath")

	# create objs and deps lists, and a list with all subdirectories to create inside OBJ_DIR
	OBJS := $(CPP_SRCS:$(SOURCE_DIR)/%.cpp=$(OBJ_DIR)/%.o)
	CPP_DEPS := $(CPP_SRCS:$(SOURCE_DIR)/%.cpp=$(OBJ_DIR)/%.d)
	TREE_WINDOWS := $(sort $(patsubst %/,%,$(dir $(OBJS))))
else
	# if on unix, search for all .cpp files from sources directory
	CPP_SRCS := $(shell find $(SOURCE_DIR) -name "*.cpp")

	# create objs and deps lists, and a list with all subdirectories to create inside OBJ_DIR
	OBJS := $(CPP_SRCS:$(SOURCE_DIR)/%.cpp=$(OBJ_DIR)/%.o)
	CPP_DEPS := $(CPP_SRCS:$(SOURCE_DIR)/%.cpp=$(OBJ_DIR)/%.d)
	TREE_UNIX := $(sort $(patsubst %/,%,$(dir $(OBJS))))
endif



# all target
all: $(BUILD_DIR)/ISSH

$(BUILD_DIR)/ISSH: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CC) -o "$@" $(OBJS) $(USER_OBJS) $(LIBS) -fopenmp
	@echo 'Finished building target: $@'
	@echo ' '


# debug target
debug: $(DEBUG_DIR)/ISSH

$(DEBUG_DIR)/ISSH: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CC) -o "$@" $(OBJS) $(USER_OBJS) $(LIBS) -fopenmp
	@echo 'Finished building target: $@'
	@echo ' '


# compile all dependencies for ISSH
.SECONDEXPANSION:
$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $$(@D)
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CC) $(CFLAGS) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


# create subdirectories tree, depending on which system we are

$(TREE_WINDOWS): %:
	MKDIR $@
	MKDIR $(@:$(OBJ_DIR)%=$(DEPS_DIR)%)

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
