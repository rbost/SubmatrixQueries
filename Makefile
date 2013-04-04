CC = gcc
CXX = g++  
LD = g++
CFLAGS = -Wall -O3
EXEC_NAME = test_queries
INCLUDES =
LIBS = -lrt
OBJ_FILES = matrix.o oracle_monge_matrix.o max_value.o range.o range_query.o envelope.o envelope_tree.o tests.o main.o

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
