cc = g++
obj = main.o ptt_readin.o score.o cJSON.o opdb.o region.o
main = ../main
x= -I/usr/include/mysql -L/usr/lib -lmysqlclient

all: $(obj)
	$(cc) -o $(main) $(obj) $(x) -o ../main

cJSON.o : cJSON/cJSON.c cJSON/cJSON.h
	$(cc) -c cJSON/cJSON.c $(x)

main.o: main.cpp main.h cJSON/cJSON.h
	$(cc) -c main.cpp $(x)

ptt_readin.o : ptt_readin.cpp main.h cJSON/cJSON.h
	$(cc) -c ptt_readin.cpp $(x)

score.o : score.cpp main.h cJSON/cJSON.h
	$(cc) -c score.cpp $(x)
	- rm ../Database/DataCache/*

opdb.o : opdb.cpp main.h cJSON/cJSON.h
	$(cc) -c opdb.cpp $(x)

region.o : region.cpp main.h cJSON/cJSON.h
	$(cc) -c region.cpp $(x)

.PHONY : clean
clean :
	rm $(obj)

.PHONY : cleanall
cleanall :
	rm $(obj) ../main

