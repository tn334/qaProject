# prompt> make
# builds everything and links in test program qaProject
# needs deepstate, with libfuzzer and afl support

SRCS = qaProject.cpp
OBJS = qaProject.o
OBJSDS = qaProject.o
OBJSDSSAN = san_qaProject.o
OBJSDSAFL = afl_qaProject.o

CC=clang++

CFLAGS = -fopenmp -lm

UNIT = my_test
#deepstate exec
DS = ds_test
#DeepState exec w sanitizers
DSSAN = ds_test_san
#AFL exec
DSAFL = ds_test_afl

all: $(UNIT) $(DS) $(DSSAN) $(DSAFL)

$(UNIT):	$(OBJS)
			$(CC) $(CFLAGS) $(OBJS) -o $(UNIT)

$(DS):		$(OBJSDS)
			$(CC) -std=c++14 $(CFLAGS) -o $(DS) deepstateTests.cpp $(OBJSDS) -ldeepstate

$(DSSAN):	$(OBJSDSSAN)
			$(CC) -std=c++14 $(CFLAGS) -fsanitize=undefined,integer,address -o $(DSSAN) deepstateTests.cpp $(OBJSDSSAN) -ldeepstate

$(DSAFL):	$(OBJSDSAFL) deepstateTests.cpp
			${AFL_HOME}/afl-clang++ -std=c++14 $(CFLAGS) -o $(DSAFL) deepstateTests.cpp $(OBJSDSAFL) -ldeepstate_LF

qaProject.o:	qaProject.cpp

san_qaProject.o:	qaProject.cpp
					$(CC) $(CFLAGS) -c -o san_qaProject.ol qaProject.cpp -fsanitize=undefined,address,integer

afl_qaProject.o:	qaProject.cpp
					${AFL_HOME}/afl-clang $(CFLAGS) -c -o afl_qaProject.o qaProject.cpp

clean:			
	rm -f *.o *~ $(UNIT) $(DS) $(DSSAN) $(DSLF) $(DSAFL) *.gcda *.gcno *.gcov
