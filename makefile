cc = g++
obj = main.o score.o cJSON.o region.o
main = ../main
x= -I/usr/include/mysql -L/usr/lib -lmysqlclient

all: $(obj)
	$(cc) -o $(main) $(obj) $(x) -o ../main

cJSON.o : cJSON/cJSON.c cJSON/cJSON.h
	$(cc) -c cJSON/cJSON.c $(x)

main.o: main.cpp main.h cJSON/cJSON.h
	$(cc) -c main.cpp $(x)

score.o : score.cpp main.h cJSON/cJSON.h
	$(cc) -c score.cpp $(x)
	- rm ../Database/DataCache/*

region.o : region.cpp main.h cJSON/cJSON.h
	$(cc) -c region.cpp $(x)

.PHONY : clean
clean :
	rm $(obj)

.PHONY : cleanall
cleanall :
	rm $(obj) ../main

