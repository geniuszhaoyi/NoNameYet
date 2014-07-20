obj = main.o ptt_readin.o score.o cJSON.o

all: $(obj)
	cc -o edit $(obj)

cJSON.o : cJSON/cJSON.c cJSON/cJSON.h
	cc -c cJSON/cJSON.cpp

main.o ptt_readin.o score.o : Main.h cJSON/cJSON.h

.PHONY : clean
clean :
	rm edit $(obj)
