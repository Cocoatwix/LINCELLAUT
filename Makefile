
#Compiling shared libraries must be done on Windows since Python is running on Windows
#Use cc -fPIC -shared -o orbitvis.so orbitvis.c to compile shared library

COMPILER = gcc
FLAGS = -Wall -Wextra -pedantic
TARGET = lincellaut

LIBPATH = libraries
LIB1 = helper
LIB2 = bigint
LIB3 = algebra
LIB4 = modular
LIB5 = factors
LIB6 = linalg
LIB7 = cycles

OBJPATH = objects
ALLOBJS = $(OBJPATH)/$(LIB1).o 
ALLOBJS += $(OBJPATH)/$(LIB2).o 
ALLOBJS += $(OBJPATH)/$(LIB3).o
ALLOBJS += $(OBJPATH)/$(LIB4).o
ALLOBJS += $(OBJPATH)/$(LIB5).o
ALLOBJS += $(OBJPATH)/$(LIB6).o
ALLOBJS += $(OBJPATH)/$(LIB7).o

all:	$(TARGET)

verbose:	FLAGS += -DVERBOSE
verbose:	$(TARGET)

#The -lm is included at the end of the command to avoid linker problems
# stackoverflow.com/questions/11336477
$(TARGET):	$(TARGET).c $(ALLOBJS)
	$(COMPILER) $(FLAGS) -o $(TARGET) $(TARGET).c $(ALLOBJS) -lm
	
$(OBJPATH)/$(LIB1).o:	$(LIBPATH)/$(LIB1).c | $(OBJPATH)
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB1).o $(LIBPATH)/$(LIB1).c
	
$(OBJPATH)/$(LIB2).o:	$(LIBPATH)/$(LIB2).c | $(OBJPATH)
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB2).o $(LIBPATH)/$(LIB2).c
	
$(OBJPATH)/$(LIB3).o:	$(LIBPATH)/$(LIB3).c | $(OBJPATH)
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB3).o $(LIBPATH)/$(LIB3).c
	
$(OBJPATH)/$(LIB4).o:	$(LIBPATH)/$(LIB4).c | $(OBJPATH)
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB4).o $(LIBPATH)/$(LIB4).c
	
$(OBJPATH)/$(LIB5).o:	$(LIBPATH)/$(LIB5).c | $(OBJPATH)
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB5).o $(LIBPATH)/$(LIB5).c
	
$(OBJPATH)/$(LIB6).o:	$(LIBPATH)/$(LIB6).c | $(OBJPATH)
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB6).o $(LIBPATH)/$(LIB6).c
	
$(OBJPATH)/$(LIB7).o:	$(LIBPATH)/$(LIB7).c | $(OBJPATH)
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB7).o $(LIBPATH)/$(LIB7).c
	
#Create objects directory if it doesn't already exist
# stackoverflow.com/questions/12605051
$(OBJPATH):
	mkdir -p $@

clean:
	rm $(TARGET)
	rm $(OBJPATH)/*.o
	rm $(OBJPATH)/*.so
	rm *.orbits
	rm *.orbitsloc
	rm *.iteration