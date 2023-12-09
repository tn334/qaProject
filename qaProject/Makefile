# prompt> make
# builds everything and links in test program qaProject
# needs deepstate, with libfuzzer and afl support

SRCS = qaProject.cpp
HDRS = qaProject.h
OBJS = qaProject.o
OBJSDS = qaProject.o
OBJSDSSAN = san_qaProject.o
OBJSDSAFL = afl_qaProject.o


ifeq ($(origin CXX),default)
CXX = clang++
endif

CFLAGS = -fopenmp
AFL_HOME = /usr/local/bin/afl-gcc

# DeepState test file
DEEPSTATE_TEST_FILE = deepstateTests.cpp

#real code
UNIT = my_code
#deepstate exec
DS = ds_test
#DeepState exec w sanitizers
DSSAN = ds_test_san
#AFL exec
DSAFL = ds_test_afl

all: $(UNIT) $(DS) $(DSSAN) $(DSAFL)

$(UNIT):	$(OBJS)
			$(CXX) $(CFLAGS) $(OBJS) -o $(UNIT) $(DMALLOC_LIB)

$(DS):		$(OBJSDS) deepstateTests.cpp
			$(CXX) -std=c++14 $(CFLAGS) -o $(DS) $(DEEPSTATE_TEST_FILE) $(OBJSDS) -ldeepstate

$(DSSAN):	$(OBJSDSSAN) deepstateTests.cpp
			$(CXX) -std=c++14 $(CFLAGS) -fsanitize=undefined,integer,address -o $(DSSAN) $(DEEPSTATE_TEST_FILE) $(OBJSDSSAN) -ldeepstate

$(DSAFL):	$(OBJSDSAFL) deepstateTests.cpp
			${AFL_HOME}/afl-clang++ -std=c++14 $(CFLAGS) -o $(DSAFL) $(DEEPSTATE_TEST_FILE) $(OBJSDSAFL) -ldeepstate_AFL

qaProject.o:	qaProject.cpp qaProject.h
				#$(CXX) -c qaProject.cpp

#deepStateTests.o: deepstateTests.cpp qaProject.cpp qaProject.h

san_qaProject.o:	qaProject.cpp qaProject.h
					$(CXX) $(CFLAGS) -c -o san_qaProject.o qaProject.cpp -fsanitize=undefined,address,integer

afl_qaProject.o:	qaProject.cpp qaProject.h
					${AFL_HOME}/afl-clang $(CFLAGS) -c -o afl_qaProject.o qaProject.cpp

clean:			
	rm -f *.o *~ $(UNIT) $(DS) $(DSSAN) $(DSLF) $(DSAFL) *.gcda *.gcno *.gcov
