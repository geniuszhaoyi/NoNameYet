cc = g++
obj = main.o ptt_readin.o score.o cJSON.o opdb.o
main = ../main

all: $(obj)
	$(cc) -o $(main) $(obj)

cJSON.o : cJSON/cJSON.c cJSON/cJSON.h
	$(cc) -c cJSON/cJSON.c

main.o: main.cpp main.h cJSON/cJSON.h
	$(cc) -c main.cpp

ptt_readin.o : ptt_readin.cpp main.h cJSON/cJSON.h
	$(cc) -c ptt_readin.cpp

score.o : score.cpp main.h cJSON/cJSON.h
	$(cc) -c score.cpp

opdb.o : opdb.cpp main.h cJSON/cJSON.h
	$(cc) -c opdb.cpp

.PHONY : clean
clean :
	rm $(obj)

.PHONY : cleanall
cleanall :
	rm $(obj) ../main

