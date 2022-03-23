all: env

env: main.o random_utilities.o Observer.o
	g++ main.o random_utilities.o Observer.o -o env

random_utilities.o: random_utilities.cpp random_utilities.h
	g++ -c random_utilities.cpp

Observer.o: Observer.cpp Observer.h
	g++ -c Observer.cpp

main.o: main.cpp Observer.h random_utilities.h
	g++ -c main.cpp

clean:
	rm -rf *.o env
