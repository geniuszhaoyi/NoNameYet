#include "main.h"

typedef struct node{
    int s,t;
    int type;
}node;

int nt;
node ns[NODE_SIZE];

bool cmp2(node a,node b){
    if(a.s==b.s) return a.type>b.type;
    return a.s<b.s;
}

int get_Chr_No(const char *specie,const char *chr_name){
    int res;
    MYSQL_RES *result;
    MYSQL_ROW sql_row;
    char buffer[200];
    sprintf(buffer,"SELECT Chr_No FROM Table_Specie as s JOIN Table_chromosome as c WHERE s.Sno=c.Sno and SName='%s' and Chr_Name='%s';",specie,chr_name);
    res=mysql_query(my_conn,buffer);
    if(res) return -1;
    result=mysql_store_result(my_conn);
    if((sql_row=mysql_fetch_row(result))){
        return atoi(sql_row[0]);
    }
    return -1;
}

int region_wmin(int a,int b){
    const int region_cw[]={0,2,3,1,4,0};
    return region_cw[a]<region_cw[b]?a:b;
}

int region_str[GENE_LEN];

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
        if(sql_row[2][0]!='N') utrstart+=atoi(sql_row[2])-1;
        if(sql_row[3][0]!='N') utrend-=atoi(sql_row[3])+1+3;
        //new node for gene
        ns[nt].s=start;
        ns[nt].t=end;
        ns[nt].type=REGION_INTRON;
        nt++;
        //new node for UTR
        if(start!=utrstart){
            ns[nt].s=start;
            ns[nt].t=utrstart;
            ns[nt].type=REGION_UTR;
            nt++;
        }
        if(end!=utrend){
            ns[nt].s=utrend;
            ns[nt].t=end;
            ns[nt].type=REGION_UTR;
            nt++;
        }
    }
    sprintf(buffer,"SELECT region_Start, region_end FROM Table_region as r JOIN Table_gene as g where r.gene_ID=g.gene_ID and Chr_No=%d and rtd_id=1 and region_Start<=%d and region_end>=%d;",Chr_No,end,start);
    res=mysql_query(my_conn,buffer);
    if(res) return NULL ;
    result=mysql_store_result(my_conn);
    while((sql_row=mysql_fetch_row(result))){
        //new node for exon
        ns[nt].s=atoi(sql_row[0]);
        ns[nt].t=atoi(sql_row[1]);
        ns[nt].type=REGION_EXON;
        nt++;
    }

    sort(ns,ns+nt,cmp2);

    int i,j;
    for(j=0;j<=end-start;j++) region_str[j]=REGION_INTERGENIC;
    for(i=0;i<nt;i++){
        for(j=ns[i].s;j<=ns[i].t;j++){
            region_str[j-start]=region_wmin(ns[i].type,region_str[j-start]);
        }
    }
    for(j=0;j<=end-start-1;j++){
        if(region_str[j]!=region_str[j+1]){
            cJSON *n=cJSON_CreateObject();
            cJSON_AddNumberToObject(n,"endpoint",j+start);
            cJSON_AddStringToObject(n,"description",region_info[region_str[j]]);
            cJSON_AddItemToArray(root,n);
        }
    }
    cJSON *n=cJSON_CreateObject();
    cJSON_AddNumberToObject(n,"endpoint",end);
    cJSON_AddStringToObject(n,"description",region_info[region_str[end-start]]);
    cJSON_AddItemToArray(root,n);

    return root;
}

int getRegion(int sgrna_ID,int Chr_No,int sgrna_start,int sgrna_end){
    char buffer[1000];
    sprintf(buffer,"SELECT rtd_id FROM Table_sgRNA WHERE sgrna_ID=%d; ",sgrna_ID);
    int res=mysql_query(my_conn,buffer);
    if(res){
    }
    MYSQL_RES *result=mysql_store_result(my_conn);
    MYSQL_ROW sql_row;
    if((sql_row=mysql_fetch_row(result))){
        if(sql_row[0]!=NULL && sql_row[0][0]!='N' && sql_row[0][0]!='n'){
            return atoi(sql_row[0]);
        }
    }

    //First: find if in gene
    sprintf(buffer,"SELECT gene_Start, gene_end, gene_UTRstart, gene_UTRend FROM Table_gene WHERE Chr_No=%d and gene_start<=%d and gene_end>=%d",Chr_No,sgrna_start,sgrna_end);
    mysql_query(my_conn,buffer);
    result=mysql_store_result(my_conn);
    if((sql_row=mysql_fetch_row(result))){
        int start=atoi(sql_row[0]),end=atoi(sql_row[1]);
        int utrstart=atoi(sql_row[0]),utrend=atoi(sql_row[1]);
        if(sql_row[2][0]!='N') utrstart+=atoi(sql_row[2])-1;
        if(sql_row[3][0]!='N') utrend-=atoi(sql_row[3])+1+3;
        if((start<=sgrna_start && utrstart>=sgrna_end) || (utrend<=sgrna_start && end>=sgrna_end)){
            int type=REGION_UTR;
            sprintf(buffer,"UPDATE Table_sgRNA SET rtd_id=%d WHERE sgrna_ID=%d",type,sgrna_ID);
            mysql_query(my_conn,buffer);
            return type;
        }
    }else{
        int type=REGION_INTERGENIC;
        sprintf(buffer,"UPDATE Table_sgRNA SET rtd_id=%d WHERE sgrna_ID=%d",type,sgrna_ID);
        mysql_query(my_conn,buffer);
        return type;
    }

    //Secondly: find if in exon
    sprintf(buffer,"SELECT region_ID FROM Table_region as r join Table_gene as g WHERE r.gene_ID=g.gene_ID and rtd_id=1 and Chr_No=%d and region_start<=%d and region_end>=%d",Chr_No,sgrna_start,sgrna_end);
    mysql_query(my_conn,buffer);
    result=mysql_store_result(my_conn);
    if((sql_row=mysql_fetch_row(result))){
        int type=REGION_EXON;
        sprintf(buffer,"UPDATE Table_sgRNA SET rtd_id=%d WHERE sgrna_ID=%d",type,sgrna_ID);
        mysql_query(my_conn,buffer);
        return type;
    }else{
        int type=REGION_INTRON;
        sprintf(buffer,"UPDATE Table_sgRNA SET rtd_id=%d WHERE sgrna_ID=%d",type,sgrna_ID);
        mysql_query(my_conn,buffer);
        return type;
    }
}
