#include "main.h"

#define REGION_GENE 1
#define REGION_UTR 2
#define REGION_EXON 3

typedef struct node{
    int s,t;
    int type;
}node;

int nt;
node ns[NODE_SIZE];

bool cmp2(node a,node b){
    if(a.s==b.s) return a.type<b.type;
    return a.s<b.s;
}

int get_Chr_No(const char *specie,const char *chr_name){
    int res;
    MYSQL_RES *result;
    MYSQL_ROW sql_row;
    char buffer[200];
    res=sprintf(buffer,"SELECT Chr_No FROM Table_Specie as s JOIN Table_chromosome as c WHERE s.Sno=c.Sno and SName='%s' and Chr_Name='%s';",specie,chr_name);
    if(res) return -1;
    result=mysql_store_result(my_conn);
    if((sql_row=mysql_fetch_row(result))){
        return atoi(sql_row[0]);
    }
    return -1;
}

cJSON *getlineregion(int Chr_No,int start,int end){

    cJSON *root=cJSON_CreateArray();
    nt=0;
    MYSQL_RES *result;
    MYSQL_ROW sql_row;
    char buffer[200];

    int res;
    sprintf(buffer,"SELECT gene_Start, gene_end, gene_UTRstart, gene_UTRend FROM Table_gene where Chr_No=%d and gene_Start<=%d and gene_end>=%d;",Chr_No,end,start);
    res=mysql_query(my_conn,buffer);
    if(res) return NULL;
    result=mysql_store_result(my_conn);
    while((sql_row=mysql_fetch_row(result))){
        int start=atoi(sql_row[0]),end=atoi(sql_row[1]);
        int utrstart=atoi(sql_row[0]),utrend=atoi(sql_row[1]);
        if(sql_row[2][0]!='N') utrstart-=atoi(sql_row[2]);
        if(sql_row[3][0]!='N') utrend+=atoi(sql_row[3]);
        //new node for gene
        ns[nt].s=start;
        ns[nt].t=end;
        ns[nt].type=REGION_GENE;
        nt++;
        //new node for UTR
        if(start!=utrstart){
            ns[nt].s=utrstart;
            ns[nt].t=start-1;
            ns[nt].type=REGION_UTR;
            nt++;
        }
        if(end!=utrend){
            ns[nt].s=end+1;
            ns[nt].t=utrend;
            ns[nt].type=REGION_UTR;
            nt++;
        }
    }
    sprintf(buffer,"SELECT gene_Start, gene_end FROM Table_region as r JOIN Table_gene as g where r.gene_ID=g.gene_ID and Chr_No=%d and region_Start<=%d and region_end>=%d;",Chr_No,start,end);
    res=mysql_query(my_conn,buffer);
    if(res) return NULL;
    result=mysql_store_result(my_conn);
    while((sql_row=mysql_fetch_row(result))){
        //new node for exon
        ns[nt].s=atoi(sql_row[0]);
        ns[nt].t=atoi(sql_row[1]);
        ns[nt].type=REGION_EXON;
        nt++;
    }

    sort(ns,ns+nt,cmp2);

    int i;
    int last=start,mark_gene=-1;
    for(i=0;i<nt;i++){
        if(last+1<ns[i].s){
            if(ns[i].s<mark_gene){
                cJSON *n=cJSON_CreateObject();
                cJSON_AddNumberToObject(n,"endpoint",ns[i].s-1);
                cJSON_AddStringToObject(n,"description","intron");
                cJSON_AddItemToArray(root,n);
            }else{
                cJSON *n=cJSON_CreateObject();
                cJSON_AddNumberToObject(n,"endpoint",ns[i].s-1);
                cJSON_AddStringToObject(n,"description","intergenic");
                cJSON_AddItemToArray(root,n);
            }
        }
        if(ns[i].type==REGION_UTR){
            cJSON *n=cJSON_CreateObject();
            cJSON_AddNumberToObject(n,"endpoint",ns[i].t);
            cJSON_AddStringToObject(n,"description","utr");
            cJSON_AddItemToArray(root,n);
            last=ns[i].t;
        }
        if(ns[i].type==REGION_EXON){
            cJSON *n=cJSON_CreateObject();
            cJSON_AddNumberToObject(n,"endpoint",ns[i].t);
            cJSON_AddStringToObject(n,"description","exon");
            cJSON_AddItemToArray(root,n);
            last=ns[i].t;
        }
        if(ns[i].type==REGION_GENE){
            mark_gene=ns[i].t;
        }
    }
    if(last<end){
        if(ns[i].s<mark_gene){
            cJSON *n=cJSON_CreateObject();
            cJSON_AddNumberToObject(n,"endpoint",ns[i].s-1);
            cJSON_AddStringToObject(n,"description","intron");
            cJSON_AddItemToArray(root,n);
        }else{
            cJSON *n=cJSON_CreateObject();
            cJSON_AddNumberToObject(n,"endpoint",ns[i].s-1);
            cJSON_AddStringToObject(n,"description","intergenic");
            cJSON_AddItemToArray(root,n);
        }
    }

    return root;
}

