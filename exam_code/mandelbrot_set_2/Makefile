
#first make file!
CC=gcc
CP= gcc -fopenmp
#folder for headers
HDIR=./include
#flags for compiling 
CFLAGS= -c -O3 -march=native -I$(HDIR) -Wall #-g -O0
PFLAGS= -c -fopenmp -O3 -march=native -I$(HDIR) -Wall #-g -O0
DFLAGS= -g -O3 -march=native -Wall
RFLAGS= -g -O3 -march=native -Wall

OBJECTSDIR= objects
#with wildcart is NOT necessary to have "./" beucase it will look for that pattern instead og the all path
SRC := $(wildcard src/*.c)
TEST:= $(wildcard test/*.c)
OBJ := $(patsubst src/%.c, $(OBJECTSDIR)/%.o,$(SRC))
OBJ_TEST := $(patsubst test/%.c, $(OBJECTSDIR)/%.o,$(TEST))

 #patsubt is a pattern substitution
release: $(OBJECTSDIR)/main.o $(OBJ)
	$(CC) $(RFLAGS) $^ -o mandelbrot.x

debug: $(OBJECTSDIR)/main.o $del(OBJ)
	$(CC) $(DFLAGS) $^ -o debug_mandelbrot.x
	brot.x

mandelbrot.x: $(OBJECTSDIR)/main.o $(OBJ)
	$(CC) $^ -o $@ #-lm 

#the "space" is TAB!the $@evaluate the name of the current target
#hello.x: $(OBJECTSDIR)/hello.o $(OBJECTSDIR)/main.o #target:dependencies(1,2,3.... n)
#	$(CC) $(OBJECTSDIR)/main.o $(OBJECTSDIR)/hello.o -o $@ #the "space" is TAB!the $@evaluate the name of the current target

$(OBJECTSDIR)/%.o: src/%.c 
	$(CC) $(CFLAGS) $^ -o $@ 

#build up the output files usually where you launch the make, but in this way you specifically asked it to refer and build up them into the objects folder!


$(OBJECTSDIR)/main.o: main.c
	$(CC) $(CFLAGS) main.c -o $@ 
#hello.o: hello.c
#	$(CC) $(CFLAGS) hello.c
clean:
	rm -rf *.o *.x objects/*.o
print:
	@echo $(SRC)
	@echo $(OBJ)

test: test_mandelbrot.x

test_mandelbrot.x: $(OBJECTSDIR)/test.o $(OBJ_TEST)
	$(CC) $^ -o $@ #-lm 

$(OBJECTSDIR)/%.o: test/%.c 
	$(CC) $(CFLAGS) $^ -o $@ 

parallel: mandelbrot_parallel.x

mandelbrot_parallel.x: $(OBJECTSDIR)/main_p.o $(OBJ)
	$(CP) -fopenmp $^ -o $@ #-lm 

$(OBJECTSDIR)/%.o: src/%.c 
	$(CP) $(PFLAGS) $^ -o $@ 

$(OBJECTSDIR)/main_p.o: main_p.c
	$(CP) $(CFLAGS) main_p.c -o $@ 