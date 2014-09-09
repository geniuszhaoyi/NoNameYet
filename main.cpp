/** @file
@brief Main function and some assistant functions.
@author Yi Zhao
*/

#include "main.h"

///Global restricts.
restrict req_restrict;

int ini;///<Total number of in_site.
site in_site[NODE_SIZE];///<Cadidate's information.

MYSQL *my_conn;///<Mysql Connection.

///Compare sites by site.score
bool cmp_in_site(site a,site b){
    return a.score>b.score;
}
///Compare sites by site.index
bool cmp_by_index(site a,site b){
    return a.index<b.index;
}

///Read charactors until reach '\\n'
int readLine(FILE *file){
    char ch;
    while(fscanf(file,"%c",&ch)==1) if(ch=='\n') return 1;
    return 0;
}

///Create an json array from array of cJSON. Source code are from cJSON MAN.
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

///Check if cadidate sgRNA's region matches request.
int check_region(int i){
    int r=in_site[i].region;
    return req_restrict.region[r];
}

///Check if cadidate sgRNA matches request RFC restriction.
int check_rfc(int i){
    char str[LEN+PAM_LEN+3];
    strcpy(str,in_site[i].nt);
    strcat(str,in_site[i].pam);
    if(req_restrict.rfc10){
        if(strstr(str,"GAATTC")) return 0;
        if(strstr(str,"TCTAGA")) return 0;
        if(strstr(str,"ACTAGT")) return 0;
        if(strstr(str,"CTGCAG")) return 0;
        if(strstr(str,"GCGGCCGC")) return 0;
    }
    if(req_restrict.rfc12){
        if(strstr(str,"GAATTC")) return 0;
        if(strstr(str,"ACTAGT")) return 0;
        if(strstr(str,"GCTAGC")) return 0;
        if(strstr(str,"CTGCAG")) return 0;
        if(strstr(str,"GCGGCCGC")) return 0;
    }
    if(req_restrict.rfc12a){
        if(strstr(str,"CAGCTG")) return 0;
        if(strstr(str,"CTCGAG")) return 0;
        if(strstr(str,"CCTAGG")) return 0;
        if(strstr(str,"TCTAGA")) return 0;
        if(strstr(str,"GCTCTTC")) return 0;
        if(strstr(str,"GAAGAGC")) return 0;
    }
    if(req_restrict.rfc21){
        if(strstr(str,"GAATTC")) return 0;
        if(strstr(str,"AGATCT")) return 0;
        if(strstr(str,"GGATCC")) return 0;
        if(strstr(str,"CTCGAG")) return 0;
    }
    if(req_restrict.rfc23){
        if(strstr(str,"GAATTC")) return 0;
        if(strstr(str,"TCTAGA")) return 0;
        if(strstr(str,"ACTAGT")) return 0;
        if(strstr(str,"CTGCAG")) return 0;
        if(strstr(str,"GCGGCCGC")) return 0;
    }
    if(req_restrict.rfc25){
        if(strstr(str,"GAATTC")) return 0;
        if(strstr(str,"TCTAGA")) return 0;
        if(strstr(str,"GCCGGC")) return 0;
        if(strstr(str,"ACCGGT")) return 0;
        if(strstr(str,"ACTAGT")) return 0;
        if(strstr(str,"CTGCAG")) return 0;
        if(strstr(str,"GCGGCCGC")) return 0;
    }
    return 1;
}

//R=A,G; M=A,C; W=A,T; S=C,G; K=G,T; Y=C,T; H=A,C,T; V=A,C,G; B=C,G,T; D=A,G,T; N=A,G,C,T
///Check if cadidate sgRNA's PAM matches request.
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

///Reverse DNA charactor. A-T & C-T.
char dna_rev_char(char ch){
    if(ch=='A' || ch=='a') return 'T';
    if(ch=='T' || ch=='t') return 'A';
    if(ch=='C' || ch=='c') return 'G';
    if(ch=='G' || ch=='g') return 'C';
    return 'N';
}

///Get the opposite DNA seqence.
char *dna_rev(char *sr,const char *s,int len){
    int i;
    for(i=0;i<len;i++){
        sr[i]=dna_rev_char(s[len-i-1]);
    }
    sr[i]=0;
    return sr;
}

///Remove all non-readable charactors from string(char*).
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

///Not remove all non-readable charactors from string(char*). @deprecated Only for debug.
char *_NomoreSpace(char *str){
    return str;
}

///Check request json.
int check_req(cJSON *request){
    int clt=3;
    int i;
    char cl[][10]={"pam","specie","rfc"};
    for(i=0;i<clt;i++){
        if(cJSON_GetObjectItem(request,cl[i])==NULL) break;
    }
    return i<clt;
}

///Print error message.
void onError(const char *msg){
    char *p;

    cJSON *root_error=cJSON_CreateObject();
    cJSON_AddNumberToObject(root_error,"status",1);
    cJSON_AddStringToObject(root_error,"message",msg);

    printf("%s\n",NomoreSpace(p=cJSON_Print(root_error)));
    free(p);
}

///Default request.
char argv_default[]="{\"specie\":\"Saccharomyces-cerevisiae\",\"length\":17,\"location\":\"NC_001144-chromosome12:1..500\",\"pam\":\"NGG\",\"rfc\":\"100010\"}";
const char *region_info[]={"","EXON","INTRON","UTR","INTERGENIC"};///<Name of region. @see REGION_EXON, REGION_INTRON, REGION_UTR, REGION_INTERGENIC, REGION_GENE

/**
@brief Root pointer for saved possible-offtarget sgRNA's information.
@see localresult.cpp
*/
localrow *localresult;

mos_pthread_mutex_t mutex;///<Mutex for thread_share_variables. @see create_thread_socre, new_thread
mos_pthread_mutex_t mutex_mysql_conn;///<Mutex for mysql_connextion(my_conn). @see my_conn
mos_sem_t sem_thread;///<Semaphore for thread. Number of thread cann't exceed MAX_SEM_THREAD. @see MAX_SEM_THREAD

/**
Main function. Include Input, output and Database create connect.
*/
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
    req_restrict.region[0]=1;
    req_restrict.region[1]=1;
    req_restrict.region[2]=1;
    req_restrict.region[3]=1;
    req_restrict.region[4]=1;

    localresult=NULL;

    mos_pthread_mutex_init(&mutex,NULL);
    mos_pthread_mutex_init(&mutex_mysql_conn,NULL);
    mos_sem_init(&sem_thread,0,MAX_SEM_THREAD);

    char *req_str=argv_default;
    if(args==2) req_str=argv[1];

    cJSON *request=cJSON_Parse(req_str);
    if(check_req(request)){
        onError("illegal args");
        return RETUEN_ERROR;
    }

    i=1;

    cJSON *cJSON_temp;
    cJSON_temp=cJSON_GetObjectItem(request,"type");
    if(cJSON_temp) req_type=cJSON_temp->valueint;
    char req_pam[PAM_LEN+1];
    strcpy(req_pam,cJSON_GetObjectItem(request,"pam")->valuestring);
    char req_specie[30];
    strcpy(req_specie,cJSON_GetObjectItem(request,"specie")->valuestring);

    double req_r1=0.65;
    cJSON_temp=cJSON_GetObjectItem(request,"r1");
    if(cJSON_temp){
        req_r1=cJSON_GetObjectItem(request,"r1")->valuedouble;
    }

    cJSON_temp=cJSON_GetObjectItem(request,"gene");
    int req_gene_start,req_gene_end;
    char req_id[100];
    if(cJSON_temp){
        onError("temporary not available");
        return RETUEN_ERROR;
    }else{
        strcpy(buffer,cJSON_GetObjectItem(request,"location")->valuestring);
        for(i=0;buffer[i]!=':' && buffer[i]!=0;i++){
            req_id[i]=buffer[i];
        }req_id[i]=0;
        sscanf(buffer+i+1,"%d..%d",&req_gene_start,&req_gene_end);
    }

    cJSON_temp=cJSON_GetObjectItem(request,"region");
    if(cJSON_temp){
        req_restrict.region[1]=(cJSON_temp->valuestring)[0]-48;
        req_restrict.region[2]=(cJSON_temp->valuestring)[1]-48;
        req_restrict.region[3]=(cJSON_temp->valuestring)[2]-48;
        req_restrict.region[4]=(cJSON_temp->valuestring)[3]-48;
    }

    char req_rfc[10];
    strcpy(req_rfc,cJSON_GetObjectItem(request,"rfc")->valuestring);
    req_restrict.rfc10=req_rfc[0]-48;
    req_restrict.rfc12=req_rfc[1]-48;
    req_restrict.rfc12a=req_rfc[2]-48;
    req_restrict.rfc21=req_rfc[3]-48;
    req_restrict.rfc23=req_rfc[4]-48;
    req_restrict.rfc25=req_rfc[5]-48;

    req_restrict.ntlength=20;
    cJSON_temp=cJSON_GetObjectItem(request,"length");
    if(cJSON_temp){
        req_restrict.ntlength=cJSON_temp->valueint;
        if(req_restrict.ntlength>20 || req_restrict.ntlength<17){
            onError("Invaild length! ");
            return RETUEN_ERROR;
        }
    }

    /*
    This part above is for read in JSON-style request.
    The result stored in req_specie, req_pam, req_gene_start, req_gene_end, rfc and so on.
    At the same time, it reads in ptt file for specie,
    and also open files.
    */

    MYSQL_ROW sql_row;
    my_conn=mysql_init(NULL);

    my_bool mb=false;
    mysql_options(my_conn,MYSQL_SECURE_AUTH,&mb);
#ifdef  _WIN32
    if(mysql_real_connect(my_conn,"127.0.0.1","root","zy19930108","CasDB",3306,NULL,0)){
#endif
#ifdef  __linux
    if(mysql_real_connect(my_conn,"127.0.0.1","igem","uestc2014!","CasDB",3306,NULL,0)){
#endif
    }else{
        sprintf(buffer,"database connect error\n$%s",mysql_error(my_conn));
        onError(buffer);
        return RETUEN_ERROR;
    }

    int res;
    sprintf(buffer,"SELECT Sno FROM Table_Specie WHERE SName=\"%s\";",req_specie);
    res=mysql_query(my_conn,buffer);
    if(res){
        onError("database select error1");
        return RETUEN_ERROR;
    }
    MYSQL_RES *result=mysql_store_result(my_conn);
    if((sql_row=mysql_fetch_row(result))){
    }else{
        onError("No such specie!");
        return RETUEN_ERROR;
    }
    mysql_free_result(result);

    sprintf(buffer,"SELECT sgrna_start, sgrna_end, sgrna_strand, sgrna_seq, sgrna_PAM, Chr_Name, sgrna_ID, Chr_No FROM view_allsgrna WHERE SName='%s' and pam_PAM='%s';",req_specie,req_pam);
    res=mysql_query(my_conn,buffer);
    if(res){
        onError("database select error2");
        return RETUEN_ERROR;
    }
    MYSQL_RES *result_t=mysql_store_result(my_conn);
    make_mysqlres_local(&localresult,result_t);
    localres_count(localresult);
    mysql_free_result(result_t);

    sprintf(buffer,"SELECT sgrna_start, sgrna_end, sgrna_strand, sgrna_seq, sgrna_PAM, Chr_Name, sgrna_ID, Chr_No FROM view_allsgrna WHERE SName='%s' and pam_PAM='%s' and Chr_Name='%s' and sgrna_start>=%d and sgrna_end<=%d;",req_specie,req_pam,req_id,req_gene_start,req_gene_end);
    res=mysql_query(my_conn,buffer);
    if(res){
        onError("database select error1");
        return RETUEN_ERROR;
    }
    result=mysql_store_result(my_conn);
    mysql_data_seek(result,0);
    while((sql_row=mysql_fetch_row(result))){
        in_site[ini].index=atoi(sql_row[0]);
        in_site[ini].strand=sql_row[2][0];
        strcpy(in_site[ini].nt,sql_row[3]);
        strcpy(in_site[ini].pam,sql_row[4]);
        in_site[ini].ot.clear();
        strcpy(in_site[ini].chromosome,sql_row[5]);
        mos_pthread_mutex_lock(&mutex_mysql_conn);
        in_site[ini].region=getRegion(atoi(sql_row[6]),atoi(sql_row[7]),atoi(sql_row[0]),atoi(sql_row[1]));
        mos_pthread_mutex_unlock(&mutex_mysql_conn);

        if(check_region(ini)==0){
            continue;
        }
        if(check_rfc(ini)==0){
            continue;
        }

        localrow lr;
        for(int i=0;i<8;i++){
            strcpy(lr.row[i],sql_row[i]);
        }

        sprintf(buffer,"SELECT result_Sspe, result_Seff, result_count, result_offtarget FROM Table_sgRNA as r JOIN Table_result as t ON r.sgrna_ID=t.sgrna_ID WHERE result_ntlength=%d and r.sgrna_ID=%s and result_offtarget IS NOT NULL",req_restrict.ntlength,lr.row[6]);
        mos_pthread_mutex_lock(&mutex_mysql_conn);
        int res=mysql_query(my_conn,buffer);
        if(res){
            printf("%s\n",mysql_error(my_conn));
        }
        MYSQL_RES *rsdc=mysql_store_result(my_conn);
        mos_pthread_mutex_unlock(&mutex_mysql_conn);
        if((sql_row=mysql_fetch_row(rsdc))){
            sscanf(sql_row[0],"%lf",&in_site[ini].Sspe_nor);
            sscanf(sql_row[1],"%lf",&in_site[ini].Seff_nor);
            in_site[ini].score=req_r1*in_site[ini].Sspe_nor+(1-req_r1)*in_site[ini].Seff_nor;
            sscanf(sql_row[2],"%d",&in_site[ini].count);
            in_site[ini].otj=cJSON_Parse(sql_row[3]);
        }else{
            mos_sem_wait(&sem_thread);
            create_thread_socre(localresult,lr,ini,req_type,req_r1);
        }

        ini++;
    }

    for(i=0;i<ini;i++){
        if(in_site[i].ntid) mos_pthread_join(in_site[i].ntid,NULL);
    }
    free_mysqlres_local(localresult);

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
    for(i=0;i<ini && i<50;i++){
        cJSON *ans=cJSON_CreateObject();
        sprintf(buffer,"#%d",i+1);
        cJSON_AddStringToObject(ans,"key",buffer);
        sprintf(buffer,"%s%s",in_site[i].nt+(LEN-req_restrict.ntlength),in_site[i].pam);
        cJSON_AddStringToObject(ans,"grna",buffer);
        sprintf(buffer,"%s:%d",in_site[i].chromosome,in_site[i].index);
        cJSON_AddStringToObject(ans,"position",buffer);
        char xs[2];
        xs[0]=in_site[i].strand;
        xs[1]=0;
        cJSON_AddStringToObject(ans,"strand",xs);
        cJSON_AddStringToObject(ans,"region",region_info[in_site[i].region]);
        cJSON_AddNumberToObject(ans,"total_score",(int)in_site[i].score);
        cJSON_AddNumberToObject(ans,"Sspe",(int)(req_r1*in_site[i].Sspe_nor));
        cJSON_AddNumberToObject(ans,"Seff",(int)((1.0-req_r1)*in_site[i].Seff_nor));
        cJSON_AddNumberToObject(ans,"count",in_site[i].count);
        cJSON_AddItemToObject(ans,"offtarget",in_site[i].otj);

        list.push_back(ans);
    }
    cJSON_AddItemToObject(root,"result",Create_array_of_anything(&(list[0]),list.size()));

#ifdef  _WIN32
    fprintf(fopen("D:/out.txt","w"),"%s\n",_NomoreSpace(argv[0]=cJSON_Print(root)));
    printf("%s\n",NomoreSpace(argv[0]));
#endif // _WIN32
#ifdef  __linux
    //printf("{\"status\":1,\"message\":\"System in maintenance\"}");
    printf("%s\n",NomoreSpace(argv[0]=cJSON_Print(root)));
#endif // __linux

    free(argv[0]);
    mysql_free_result(result);
    mysql_close(my_conn);

    return RETURN_SUCCEED;
}
