#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
using namespace std;

#include "cJSON/cJSON.h"

#define PTT_SARS 0
#define PTT_ECOLI 1

#define LEN 20
#define PAM_LEN 20
#define NUM_NO 4

#ifdef EXPORT_HELLO_DLL
#define HELLO_API __declspec(dllexport)
#else
#define HELLO_API __declspec(dllimport)
#endif

typedef struct ptt{
    int s,t;
    char strand;
    char gene[20];
}ptt;

typedef struct site{
    char nt[LEN+1];
    char pam[PAM_LEN+1];
    int index;
    int count;
    char region;
    char strand;
    double score;
    vector <int> ot;
}site;

typedef struct restrict{
    char rfc10;
    char rfc12,rfc12a;
    char rfc21;
    char rfc23;
    char rfc25;
}restrict;

int ptt_readin(int, ptt*);

extern "C"{
    HELLO_API char *test(char *,int);
}

int readLine(FILE *);

double score(int ii,int *pini);
