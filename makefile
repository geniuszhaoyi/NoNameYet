cc = g++
obj = main.o score.o cJSON.o region.o localresult.o
main = ../main
x= -I/usr/include/mysql -L/usr/lib -lmysqlclient -pthread

all: $(obj)
	$(cc) -o $(main) $(obj) $(x)

cJSON.o : cJSON/cJSON.c cJSON/cJSON.h
	$(cc) -c cJSON/cJSON.c $(x)

main.o: main.cpp main.h cJSON/cJSON.h
	$(cc) -c main.cpp $(x)

score.o : score.cpp main.h cJSON/cJSON.h
	$(cc) -c score.cpp $(x)

region.o : region.cpp main.h cJSON/cJSON.h
	$(cc) -c region.cpp $(x)

localresult.o : localresult.cpp main.h cJSON/cJSON.h
	$(cc) -c localresult.cpp $(x)

.PHONY : clean
clean :
	rm $(obj)

.PHONY : cleanall
cleanall :
	rm $(obj) ../main

