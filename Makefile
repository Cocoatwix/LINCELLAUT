
COMPILER = gcc
FLAGS = -Wall -Wextra -pedantic
TARGET = lincellaut

LIBPATH = libraries
LIB1 = linalg
LIB2 = cycles

OBJPATH = objects
ALLOBJS = $(OBJPATH)/$(LIB1).o $(OBJPATH)/$(LIB2).o

all:	$(TARGET)

$(TARGET):	$(TARGET).c $(ALLOBJS)
	$(COMPILER) $(FLAGS) -o $(TARGET) $(TARGET).c $(ALLOBJS)
	
$(OBJPATH)/$(LIB1).o:	$(LIBPATH)/$(LIB1).c
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB1).o $(LIBPATH)/$(LIB1).c
	
$(OBJPATH)/$(LIB2).o:	$(LIBPATH)/$(LIB2).c
	$(COMPILER) $(FLAGS) -c -o $(OBJPATH)/$(LIB2).o $(LIBPATH)/$(LIB2).c
	
clean:
	rm $(OBJPATH)/*.o
	rm $(TARGET)