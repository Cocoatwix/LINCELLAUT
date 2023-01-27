
#Compiling shared libraries must be done on Windows since Python is running on Windows
#Use cc -fPIC -shared -o orbitvis.so orbitvis.c to compile shared library

COMPILER = gcc
FLAGS = -Wall -Wextra -pedantic
TARGET = lincellaut

LIBPATH = libraries
OBJPATH = objects

LIBRARIES = helper bigint algebra modular factors linalg cycles

#makefiletutorial.com/#string-substitution
FULLPATHOBJECTS = $(patsubst %,$(OBJPATH)/%.o,$(LIBRARIES))
#FULLPATHOBJECTS := $(LIBRARIES:%=$(OBJPATH)/%.o)


all:	$(TARGET)

verbose:	FLAGS += -DVERBOSE
verbose:	$(TARGET)

#valgrind --leak-check=yes --track-origins=yes -s ./lincellaut
memdebug:	FLAGS += -g -O0
memdebug:	$(TARGET)

#The -lm is included at the end of the command to avoid linker problems
# stackoverflow.com/questions/11336477
$(TARGET):	$(TARGET).c $(FULLPATHOBJECTS)
	$(COMPILER) $(FLAGS) -o $(TARGET) $(TARGET).c $(FULLPATHOBJECTS) -lm

#Generic rule for compiling all required objects
# makefiletutorial.com/#pattern-rules
# $< is the first prerequisite
# $@ is the target
$(OBJPATH)/%.o	:	$(LIBPATH)/%.c | $(OBJPATH)
	$(COMPILER) $(FLAGS) -c -o $@ $<
	
#Create objects directory if it doesn't already exist
# stackoverflow.com/questions/12605051
$(OBJPATH):
	mkdir -p $@

#makefiletutorial.com/#check-if-a-variable-is-empty
clean:
#Remove program executable if it exists
ifneq ($(TARGET),)
	rm $(TARGET)
endif
#Remove object files if they exist
ifneq ($(wildcard $(OBJPATH)/*.o),)
	rm $(wildcard $(OBJPATH)/*.o)
endif
#Remove shared object files (for ORBITVIS) if they exist
ifneq ($(wildcard $(OBJPATH)/*.so),)
	rm $(wildcard $(OBJPATH)/*.so)
endif
