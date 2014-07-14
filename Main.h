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
#define NUM_NO 7

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

int ptt_readin(int, ptt*);

int readLine(FILE *);
