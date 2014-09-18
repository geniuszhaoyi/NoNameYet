/** @file
@brief Head file for whole system.
@author Yi Zhao
@note
Include all necessary files, most definitions and structure definitions, as well as global variables.
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <map>
using namespace std;

#include "cJSON/cJSON.h"
#include "mysql.h"

#include "util.h"

#define PTT_SARS 0
#define PTT_ECOLI 1
#define PTT_SACCHAROMYCETES 2

#define REGION_EXON 1
#define REGION_INTRON 2
#define REGION_UTR 3
#define REGION_INTERGENIC 4
#define REGION_GENE 5

#define LEN 20
#define PAM_LEN 20
#define NUM_NO 4
#define NUM_CHROMOSOME 30
#define GENE_LEN 8000000
#define DCFILE_LEN 8000000
#define NODE_SIZE 1000000
#define LOCALROW_LEN 120
#define MAX_SEM_THREAD 80

#define RETUEN_ERROR -1
#define RETURN_SUCCEED 0

typedef struct site{
    char nt[LEN+1];
    char pam[PAM_LEN+1];
    int index; // position
    int count;
    char chromosome[100];
    int region;
    char strand;
    double gc;
    double score,Sspe_nor,Seff_nor;
    vector <cJSON*> ot;
    cJSON *otj;
    mos_pthread_t ntid;
}site;

typedef struct localrow{
// sgrna_start, sgrna_end, sgrna_strand, sgrna_seq, sgrna_PAM, Chr_Name, sgrna_ID, Chr_No
    char row[8][LOCALROW_LEN];
    struct localrow *next;
}localrow;

typedef struct restrict{
    char rfc10;
    char rfc12,rfc12a;
    char rfc21;
    char rfc23;
    char rfc25;
    char region[5+1];
    int ntlength;
}restrict;

struct return_struct{
    int ptts_num;
    int len[NUM_CHROMOSOME];
    int num_chromosome;
    double dou[3];
};

extern mos_pthread_mutex_t mutex;
extern mos_pthread_mutex_t mutex_mysql_conn;
extern mos_sem_t sem_thread;

extern restrict req_restrict;

extern int ini;
extern site in_site[NODE_SIZE];
extern const char *region_info[];

extern cJSON *dc_root;

extern MYSQL *my_conn;

int readLine(FILE *);
cJSON *Create_array_of_anything(cJSON **objects,int num);

void create_thread_socre(localrow *lr,localrow row,int ini,int type,double r1);

char *NomoreSpace(char *str);
char *_NomoreSpace(char *str);

int get_gene_info(char*,const char*,const char*);

int get_Chr_No(const char*,const char*);
cJSON *getlineregion(int,int,int);
int getRegion(int sgrna_ID,int Chr_No,int sgrna_start,int sgrna_end);

int make_mysqlres_local(localrow **localresult,MYSQL_RES *result_t);
void free_mysqlres_local(localrow *localresult);
int localres_count(localrow *lr);
