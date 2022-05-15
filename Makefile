all: a.out

a.out: test.c delib.c delib.h
	gcc test.c delib.c -o a.out -lm

clean:
	rm -rf *.o a.out
