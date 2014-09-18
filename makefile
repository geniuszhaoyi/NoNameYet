cc = g++
main = ../main
obj = main.o score.o cJSON.o region.o gene.o localresult.o util.o
x = -I/home/hj/immunet/www/iGEM2014/mysql-connector/include
y = -L/home/hj/immunet/www/iGEM2014/mysql-connector/lib -Wl,-dn -lmysqlclient -Wl,-dy -lm -lz -lcrypt -lpthread -ldl -lrt
util = util.h main.h cJSON/cJSON.h

all: $(obj)
	$(cc) -o $(main) $(obj) $(y)

cJSON.o : cJSON/cJSON.c cJSON/cJSON.h
	$(cc) -c cJSON/cJSON.c $(x)

util.o : util.cpp $(util)
	$(cc) -c util.cpp $(x)

main.o: main.cpp $(util)
	$(cc) -c main.cpp $(x)

score.o : score.cpp $(util)
	$(cc) -c score.cpp $(x)

region.o : region.cpp $(util)
	$(cc) -c region.cpp $(x)

gene.o : gene.cpp $(util)
	$(cc) -c gene.cpp $(x)

localresult.o : localresult.cpp $(util)
	$(cc) -c localresult.cpp $(x)

.PHONY : clean
clean :
	rm $(obj)

.PHONY : cleanall
cleanall :
	rm $(obj) ../main

