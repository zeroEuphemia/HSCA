objs = Options.o SpecificationFile.o ConstraintFile.o SAT.o main.o
CC = g++
CFLAGS = --std=c++0x -DNDEBUG -O2

all : FormatConverter

FormatConverter : $(objs)
	$(CC) $(CFLAGS) -o $@ $^

main.o : main.cc
	$(CC) $(CFLAGS) -c -o $@ $^

SAT.o : SAT.C
	$(CC) $(CFLAGS) -c -o $@ $^

ConstraintFile.o : ConstraintFile.C
	$(CC) $(CFLAGS) -c -o $@ $^

Options.o : Options.cc
	$(CC) $(CFLAGS) -c -o $@ $^

SpecificationFile.o : SpecificationFile.cc
	$(CC) $(CFLAGS) -c -o $@ $^

clean :
	-rm *.o
	-rm FormatConverter
