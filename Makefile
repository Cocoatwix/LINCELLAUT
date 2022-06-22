
#Compiling shared libraries must be done on Windows since Python is running on Windows
#Use cc -fPIC -shared -o orbitvis.so orbitvis.c to compile shared library

COMPILER = gcc
FLAGS = -Wall -Wextra -pedantic
TARGET = lincellaut

LIBPATH = libraries
LIB1 = linalg
LIB2 = cycles
LIB3 = factors
LIB4 = modular
LIB5 = fibonacci
LIB6 = bigint
LIB7 = algebra

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
# https://stackoverflow.com/questions/11336477
$(TARGET):	$(TARGET).c $(ALLOBJS)
	$(COMPILER) $(FLAGS) -o $(TARGET) $(TARGET).c $(ALLOBJS) -lm
	
$(OBJPATH)/$(LIB1).o:	$(LIBPATH)/$(LIB1).c
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB1).o $(LIBPATH)/$(LIB1).c
	
$(OBJPATH)/$(LIB2).o:	$(LIBPATH)/$(LIB2).c
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB2).o $(LIBPATH)/$(LIB2).c
	
$(OBJPATH)/$(LIB3).o:	$(LIBPATH)/$(LIB3).c
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB3).o $(LIBPATH)/$(LIB3).c
	
$(OBJPATH)/$(LIB4).o:	$(LIBPATH)/$(LIB4).c
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB4).o $(LIBPATH)/$(LIB4).c
	
$(OBJPATH)/$(LIB5).o:	$(LIBPATH)/$(LIB5).c
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB5).o $(LIBPATH)/$(LIB5).c
	
$(OBJPATH)/$(LIB6).o:	$(LIBPATH)/$(LIB6).c
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB6).o $(LIBPATH)/$(LIB6).c
	
$(OBJPATH)/$(LIB7).o:	$(LIBPATH)/$(LIB7).c
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB7).o $(LIBPATH)/$(LIB7).c

	
clean:
	rm $(TARGET)
	rm $(OBJPATH)/*.o
	rm $(OBJPATH)/*.so
	rm *.orbits
	rm *.orbitsloc
	rm *.iteration