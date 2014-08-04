#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <map>
using namespace std;

#include "cJSON/cJSON.h"
#include "mysql.h"

#define PTT_SARS 0
#define PTT_ECOLI 1
#define PTT_SACCHAROMYCETES 2

#define REGION_GENE 5
#define REGION_UTR 3
#define REGION_EXON 1
#define REGION_INTRON 2
#define REGION_INTERGENIC 4

#define LEN 20
#define PAM_LEN 20
#define NUM_NO 4
#define NUM_CHROMOSOME 30
#define GENE_LEN 8000000
#define DCFILE_LEN 8000000
#define NODE_SIZE 1000000

typedef struct site{
    char nt[LEN+1];
    char pam[PAM_LEN+1];
    int index; // position
    int count;
    char chromosome[100];
    char region;
    char strand;
    double score,Sspe_nor,Seff_nor;
    vector <cJSON*> ot;
    cJSON *otj;
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
    double dou[3];
};

extern restrict req_restrict;

extern int ini;
extern site in_site[NODE_SIZE];

extern cJSON *dc_root;

extern MYSQL *my_conn;

int readLine(FILE *);
cJSON *Create_array_of_anything(cJSON **objects,int num);

double subscore(int,int,int*,int);
void score(MYSQL_RES *result,MYSQL_ROW row,int *pini,int type,double r1);

char *NomoreSpace(char *str);
char *_NomoreSpace(char *str);

int get_Chr_No(const char*,const char*);
cJSON *getlineregion(int,int,int);
