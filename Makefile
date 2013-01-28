CC = gcc
CXX = g++  
LD = g++
CFLAGS = -Wall
EXEC_NAME = test_queries
INCLUDES =
LIBS =
OBJ_FILES = matrix.o max_value.o envelope.o envelope_tree.o main.o

all : $(EXEC_NAME)

clean :
	rm $(EXEC_NAME) $(OBJ_FILES)

$(EXEC_NAME) : $(OBJ_FILES)
	# $(CXX) -o $(EXEC_NAME) $(OBJ_FILES) $(LIBS)
	$(LD) $(OBJ_FILES) $(LIBS) -o $(EXEC_NAME)

%.o: %.cpp %.h
	$(CXX) $(CFLAGS) $(INCLUDES) -o $@ -c $<

%.o: %.cc %.h
	$(CXX) $(CFLAGS) $(INCLUDES) -o $@ -c $<

%.o: %.c %.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
