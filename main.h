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

#define LEN 20
#define PAM_LEN 20
#define NUM_NO 4
#define NUM_CHROMOSOME 30
#define GENE_LEN 8000000
#define DCFILE_LEN 8000000
#define NODE_SIZE 1000000

typedef struct ptt{
    int s,t;
    int chromosome;
    char strand;
    char gene[20];
}ptt;

typedef struct site{
    char nt[LEN+1];
    char pam[PAM_LEN+1];
    int index; // position
    int count;
    int chromosome;
    char region;
    char strand;
    double score,Sspe_nor,Seff_nor;
    vector <int> ot;
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

extern ptt ptts[NODE_SIZE];

extern int pi;
extern site psb_site[NODE_SIZE];

extern int ini;
extern site in_site[NODE_SIZE];

extern char str[NUM_CHROMOSOME][GENE_LEN];

extern cJSON *dc_root;

extern MYSQL *my_conn;

struct return_struct info_readin(int type,ptt *ptts,char str[][GENE_LEN]);

int readLine(FILE *);
cJSON *Create_array_of_anything(cJSON **objects,int num);

double subscore(int,int,int*,int);
return_struct score(int,int*,int,double);

char *NomoreSpace(char *str);
char *_NomoreSpace(char *str);

void generate_filename(char *str,const char *req_specie,const char *pam,const int type);
void dc_init(const char *fn);
void dc_save();
cJSON *dc_put(int islegal,int ini);
cJSON *dc_get(int gene,int position,char strand);

int get_Chr_No(const char*,const char*);
cJSON *getlineregion(int,int,int);
