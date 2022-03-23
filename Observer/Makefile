all: env

env: main.o Observer.o
	g++ main.o Observer.o -o env

Observer.o: Observer.cpp Observer.h
	g++ -c Observer.cpp

main.o: main.cpp Observer.h
	g++ -c main.cpp

clean:
	rm -rf *.o env