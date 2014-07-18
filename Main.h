#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
using namespace std;

#include "cJSON/cJSON.h"

#define PTT_SARS 0
#define PTT_ECOLI 1
#define PTT_SACCHAROMYCETES 2

#define LEN 20
#define PAM_LEN 20
#define NUM_NO 4
#define NUM_CHROMOSOME 30
#define GENE_LEN 8000000

typedef struct ptt{
    int s,t;
    int chromosome;
    char strand;
    char gene[20];
}ptt;

typedef struct site{
    char nt[LEN+1];
    char pam[PAM_LEN+1];
    int index;
    int count;
    int chromosome;
    char region;
    char strand;
    int score;
    vector <int> ot;
}site;

typedef struct restrict{
    char rfc10;
    char rfc12,rfc12a;
    char rfc21;
    char rfc23;
    char rfc25;
}restrict;

struct return_struct{
    int ptts_num;
    int len[NUM_CHROMOSOME];
    int num_chromosome;
};

extern restrict req_restrict;

extern ptt ptts[1000000];

extern int pi;
extern site psb_site[1000000];

extern int ini;
extern site in_site[1000000];

extern char str[NUM_CHROMOSOME][GENE_LEN];
extern char wai[NUM_CHROMOSOME][GENE_LEN];

struct return_struct info_readin(int,ptt*,char[][GENE_LEN],char[][GENE_LEN],const char*);

int readLine(FILE *);

double subscore(int,int,int*,int);
double score(int,int*);
