#include "main.h"

//req_GobalParameter
restrict req_restrict;

int ini;
site in_site[1000000];

MYSQL *my_conn;

bool cmp_in_site(site a,site b){
    return a.score>b.score;
}
bool cmp_by_index(site a,site b){
    return a.index<b.index;
}

int readLine(FILE *file){
    char ch;
    while(fscanf(file,"%c",&ch)==1) if(ch=='\n') return 1;
    return 0;
}

cJSON *Create_array_of_anything(cJSON **objects,int num)
{
	int i;cJSON *prev, *root=cJSON_CreateArray();
	for (i=0;i<num;i++)
	{
		if (!i)	root->child=objects[i];
		else	prev->next=objects[i], objects[i]->prev=prev;
		prev=objects[i];
	}
	return root;
}

//R=A,G; M=A,C; W=A,T; S=C,G; K=G,T; Y=C,T; H=A,C,T; V=A,C,G; B=C,G,T; D=A,G,T; N=A,G,C,T
int check_pam(const char *str,const char *pam){
    for(;*pam;pam++,str++){
        if(*pam=='R' && (*str=='A' || *str=='G')) continue;
        if(*pam=='M' && (*str=='A' || *str=='C')) continue;
        if(*pam=='W' && (*str=='A' || *str=='T')) continue;
        if(*pam=='S' && (*str=='C' || *str=='G')) continue;
        if(*pam=='K' && (*str=='G' || *str=='T')) continue;
        if(*pam=='Y' && (*str=='C' || *str=='T')) continue;
        if(*pam=='H' && (*str=='A' || *str=='C' || *str=='T')) continue;
        if(*pam=='V' && (*str=='A' || *str=='C' || *str=='G')) continue;
        if(*pam=='B' && (*str=='C' || *str=='G' || *str=='T')) continue;
        if(*pam=='D' && (*str=='A' || *str=='G' || *str=='T')) continue;
        if(*pam=='N' || *pam=='n') continue;
        if(*pam==*str) continue;
        return 0;
    }
    return 1;
}

char dna_rev_char(char ch){
    if(ch=='A' || ch=='a') return 'T';
    if(ch=='T' || ch=='t') return 'A';
    if(ch=='C' || ch=='c') return 'G';
    if(ch=='G' || ch=='g') return 'C';
    return 'N';
}

char *dna_rev(char *sr,const char *s,int len){
    int i;
    for(i=0;i<len;i++){
        sr[i]=dna_rev_char(s[len-i-1]);
    }
    sr[i]=0;
    return sr;
}

char *NomoreSpace(char *str){
    int i,j;
    int mark=0;
    for(i=0,j=0;str[i];i++){
        if(str[i]=='\"') mark=1-mark;
        if(mark || (str[i]!=' ' && str[i]!='\n' && str[i]!='\t')) str[j++]=str[i];
    }
    str[j]=0;
    return str;
}

char *_NomoreSpace(char *str){
    return str;
}

int check_req(cJSON *request){
    int clt=3;
    int i;
    char cl[][10]={"pam","specie","rfc"};
    for(i=0;i<clt;i++){
        if(cJSON_GetObjectItem(request,cl[i])==NULL) break;
    }
    return i<clt;
}

void onError(const char *msg){
    char *p;

    cJSON *root_error=cJSON_CreateObject();
    cJSON_AddNumberToObject(root_error,"status",1);
    cJSON_AddStringToObject(root_error,"message",msg);

    printf("%s\n",NomoreSpace(p=cJSON_Print(root_error)));
    free(p);
}

char argv_default[]="{\"specie\":\"Saccharomyces cerevisiae\",\"location\":\"NC_001134-chromosome2:200..2873\",\"pam\":\"NGG\",\"rfc\":\"100010\"}";

int main(int args,char *argv[]){
    int i;
    cJSON *root;
    int req_type=1;
    char buffer[9182+1];
    char req_gene[20]={0};

    ini=0;
    req_restrict.rfc10=0;
    req_restrict.rfc12=0;
    req_restrict.rfc12a=0;
    req_restrict.rfc21=0;
    req_restrict.rfc23=0;
    req_restrict.rfc25=0;

    char *req_str=argv_default;
    if(args==2) req_str=argv[1];

    cJSON *request=cJSON_Parse(req_str);
    if(check_req(request)){
        onError("illegal args");
        return 0;
    }

    i=1;

    cJSON *cJSON_temp;
    cJSON_temp=cJSON_GetObjectItem(request,"type");
    if(cJSON_temp) req_type=cJSON_temp->valueint;
    char req_pam[PAM_LEN+1];
    strcpy(req_pam,cJSON_GetObjectItem(request,"pam")->valuestring);
    char req_specie[30];
    strcpy(req_specie,cJSON_GetObjectItem(request,"specie")->valuestring);
    if(strcmp(req_specie,"Saccharomyces cerevisiae")==0){
    }else{
        onError("no specie");
        return 0;
    }

    double req_r1=0.65;
    cJSON_temp=cJSON_GetObjectItem(request,"gene");
    if(cJSON_temp){
        req_r1=cJSON_GetObjectItem(request,"r1")->valuedouble;
    }

    cJSON_temp=cJSON_GetObjectItem(request,"gene");
    int req_gene_start,req_gene_end;
    char req_id[100];
    if(cJSON_temp){
        onError("temporary not available");
        return 0;
    }else{
        strcpy(buffer,cJSON_GetObjectItem(request,"location")->valuestring);
        for(i=0;buffer[i]!=':' && buffer[i]!=0;i++){
            req_id[i]=buffer[i];
        }req_id[i]=0;
        sscanf(buffer+i+1,"%d..%d",&req_gene_start,&req_gene_end);
    }

    char req_rfc[10];
    strcpy(req_rfc,cJSON_GetObjectItem(request,"rfc")->valuestring);
    req_restrict.rfc10=req_rfc[0]-48;
    req_restrict.rfc12=req_rfc[1]-48;
    req_restrict.rfc12a=req_rfc[2]-48;
    req_restrict.rfc21=req_rfc[3]-48;
    req_restrict.rfc23=req_rfc[4]-48;
    req_restrict.rfc25=req_rfc[5]-48;

    /*
    This part above is for read in JSON-style request.
    The result stored in req_specie, req_pam, req_gene_start, req_gene_end, rfc and so on.
    At the same time, it reads in ptt file for specie,
    and also open files.
    */

    MYSQL_ROW sql_row;
    my_conn=mysql_init(NULL);
    if(mysql_real_connect(my_conn,"127.0.0.1","root","zy19930108","db",3306,NULL,0)){
    }else{
        sprintf(buffer,"database connect error\n$%s",mysql_error(my_conn));
        onError(buffer);
        return 0;
    }

    int res;
    sprintf(buffer,"SELECT sgrna_start, sgrna_end, sgrna_strand, sgrna_seq, sgrna_PAM, Chr_Name FROM view_getsgrna WHERE SName='%s' and sgrna_PAM='%s' and Chr_Name='%s' and sgrna_start>=%d and sgrna_end<=%d;",req_specie,req_pam,req_id,req_gene_start,req_gene_end);
    res=mysql_query(my_conn,buffer);
    if(res){
        onError("database select error1");
        return 0;
    }
    MYSQL_RES *result=mysql_store_result(my_conn);
    sprintf(buffer,"SELECT sgrna_start, sgrna_end, sgrna_strand, sgrna_seq, sgrna_PAM, Chr_Name FROM view_getsgrna WHERE SName='%s' and sgrna_PAM='%s';",req_specie,req_pam);
    res=mysql_query(my_conn,buffer);
    if(res){
        onError("database select error2");
        return 0;
    }
    MYSQL_RES *result_t=mysql_store_result(my_conn);

    mysql_data_seek(result,0);
    while((sql_row=mysql_fetch_row(result))){
        score(result_t,sql_row,&ini,req_type,req_r1);
    }

    sort(in_site,in_site+ini,cmp_in_site);  // Sort & Output

    root=cJSON_CreateObject();
    cJSON_AddNumberToObject(root,"status",0);

    cJSON_AddStringToObject(root,"specie",req_specie);
    cJSON_AddStringToObject(root,"gene",req_gene);
    sprintf(buffer,"%s:%d..%d",req_id,req_gene_start,req_gene_end);
    cJSON_AddStringToObject(root,"location",buffer);
    cJSON_temp=getlineregion(get_Chr_No(req_specie,req_id),req_gene_start,req_gene_end);    //temporary change
    if(cJSON_temp) cJSON_AddItemToObject(root,"region",cJSON_temp);

    vector<cJSON*> list;
    list.clear();
    for(i=0;i<ini;i++){
        cJSON *ans=cJSON_CreateObject();
        sprintf(buffer,"#%d",i+1);
        cJSON_AddStringToObject(ans,"key",buffer);
        sprintf(buffer,"%s%s",in_site[i].nt,in_site[i].pam);
        cJSON_AddStringToObject(ans,"grna",buffer);
        sprintf(buffer,"%s:%d",in_site[i].chromosome,in_site[i].index);
        cJSON_AddStringToObject(ans,"position",buffer);
        char xs[2];
        xs[0]=in_site[i].strand;
        xs[1]=0;
        cJSON_AddStringToObject(ans,"strand",xs);
        cJSON_AddNumberToObject(ans,"total_score",(int)in_site[i].score);
        cJSON_AddNumberToObject(ans,"Sspe",(int)(req_r1*in_site[i].Sspe_nor));
        cJSON_AddNumberToObject(ans,"Seff",(int)((1.0-req_r1)*in_site[i].Seff_nor));
        cJSON_AddNumberToObject(ans,"count",in_site[i].count);
        cJSON_AddItemToObject(ans,"offtarget",in_site[i].otj);

        list.push_back(ans);
    }
    cJSON_AddItemToObject(root,"result",Create_array_of_anything(&(list[0]),list.size()));

    fprintf(fopen("D:/out.txt","w"),"%s\n",_NomoreSpace(argv[0]=cJSON_Print(root)));
    //printf("%s\n",_NomoreSpace(argv[0]=cJSON_Print(root)));

    free(argv[0]);
    mysql_free_result(result);
    mysql_close(my_conn);

    return 0;
}
