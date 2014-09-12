/** @file
@brief Include all functions about scoring and counting.
@author Yi Zhao
*/
#include "main.h"

///Weight Values. Reference: Genetic Screens in Human Cells Using the CRISPR-Cas9 System,Wang et al., 2014
const double M[]={0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583};
///Weight Values e^(1-M[n]).
const double eM[]={2.718281828,2.718281828,2.680491036,2.718281828,2.718281828,1.831252209,1.979808257,2.718281828,1.842272751,2.511800935,1.741940985,1.635584119,1.472556491,1.160672989,1.30734714,1.187677833,1.469614321,1.216526905,1.370259311,1.517402513};

///Compare two sgRNA's Smm in reverse order.
bool cmp(cJSON *a,cJSON *b){
    double as=cJSON_GetObjectItem(a,"oscore")->valuedouble;
    double bs=cJSON_GetObjectItem(b,"oscore")->valuedouble;
    return as>bs;
}

//MYSQL_ROW row :  0:sgrna_start, 1:sgrna_end, 2:sgrna_strand, 3:sgrna_seq, 4:sgrna_PAM, 5:Chr_Name, 6:sgrna_ID, 7:Chr_No
/**
@brief Package possible-offtarget sgRNA's Infomation to cJSON Object.
@param lr possible-offtarget sgRNA's infomation saved in structure localrow.
@param oscore sgRNA's score(Smm).
@param omms Number of mismatches.
@see localrow
*/
cJSON *cJSON_otj(localrow *lr,double oscore,int omms){
    char buffer[4096];
    cJSON *root=cJSON_CreateObject();
    sprintf(buffer,"%s%s",lr->row[3]+(LEN-req_restrict.ntlength),lr->row[4]);
    cJSON_AddStringToObject(root,"osequence",buffer);
    sprintf(buffer,"%.2f",oscore);
    cJSON_AddStringToObject(root,"oscore",buffer);
    cJSON_AddNumberToObject(root,"omms",omms);
    sprintf(buffer,"%s:%s",lr->row[5],lr->row[0]);
    cJSON_AddStringToObject(root,"oposition",buffer);
    cJSON_AddStringToObject(root,"ostrand",lr->row[2]);
    mos_pthread_mutex_lock(&mutex_mysql_conn);
    cJSON_AddStringToObject(root,"oregion",region_info[getRegion(atoi(lr->row[6]),atoi(lr->row[7]),atoi(lr->row[0]),atoi(lr->row[1]))]);
    mos_pthread_mutex_unlock(&mutex_mysql_conn);
    return root;
}
/**
@brief Calculate score of a possible-offtarget sgRNA (Smm).
@param ini ID of the candidate sgRNA.
@param lr possible-offtarget sgRNA's infomation saved in structure localrow.
@param Nph Nph(type=1) or nmm(type=0).
@param type Nph's type.
@see localrow, site, in_site
*/
double subscore(int ini,localrow *lr,int *Nph,int type){
    int nmm=0;
    int d0=0;
    double smm=0;
    int nph=0;
    int i;
    for(i=LEN-req_restrict.ntlength;i<LEN;i++){
        if(in_site[ini].nt[i]!=lr->row[3][i]){
            smm+=eM[i];
            d0+=i-(LEN-req_restrict.ntlength);
            nmm++;
        }
    }
    if(nmm==0){
        nph++;
        smm=25.0;
        if(type==1) in_site[ini].ot.push_back(cJSON_otj(lr,smm,nmm));
    }else{
        if(nmm<=NUM_NO){
            smm=4.0*smm/(double)nmm/(double)nmm/(4.0*d0/((double)req_restrict.ntlength-1.0)/(double)nmm+1);
            if(type==1) in_site[ini].ot.push_back(cJSON_otj(lr,smm,nmm));
        }else{
            smm=0.0;
        }
    }
    if(type==0) (*Nph)=nmm;
    if(type==1) (*Nph)+=nph;
    return smm;
}

/**
@brief Calculate score of a candidate sgRNA.
@param lr root pointer of saved possible-offtarget sgRNA in localrow linked list.
@param row node of candidate sgRNA's infomation.
@param ini ID of the candidate sgRNA.
@param type Calculate type.
@param r1 Weight of Sspe.
@see localrow, site, in_site
*/
void score(localrow *lr,localrow row,int ini,int type,double r1){
    double r2=1.0-r1;
    int Sgc=0,S20=0;
    int Nph=0;
    double sum=0;
    int gc=0;
    int i;
    char buffer[9182];

    mos_pthread_mutex_unlock(&mutex);
    for(i=LEN-req_restrict.ntlength;i<LEN;i++) if(in_site[ini].nt[i]=='C' || in_site[ini].nt[i]=='G') gc++;
    in_site[ini].gc=(double)gc/(double)req_restrict.ntlength;
    if(in_site[ini].gc<0.4 || in_site[ini].gc>0.8) Sgc=65;
    else if(in_site[ini].gc>0.5 && in_site[ini].gc<0.7) Sgc=0;
    else Sgc=35;
    if(in_site[ini].nt[LEN-req_restrict.ntlength]!='G') S20=35;

    while(lr){
        int start=atoi(lr->row[0]);
        if(in_site[ini].index!=start){
            double smm=subscore(ini,lr,&Nph,1);
            sum+=smm;
        }
        lr=lr->next;
    }
    //sum=sigma+S1
    if(type==1 && Nph>3){
        in_site[ini].Sspe_nor=max(100-sum,0.0);
        in_site[ini].Seff_nor=100-Sgc-S20;
        in_site[ini].score=r1*in_site[ini].Sspe_nor+r2*in_site[ini].Seff_nor;
        in_site[ini].count=in_site[ini].ot.size();
    }else if(type==1){
        in_site[ini].Sspe_nor=max(100-sum,0.0);
        in_site[ini].Seff_nor=100-Sgc-S20;
        in_site[ini].score=r1*in_site[ini].Sspe_nor+r2*in_site[ini].Seff_nor;
        in_site[ini].count=in_site[ini].ot.size();
    }else{
        sum=sum-Sgc-S20+7;
        in_site[ini].Sspe_nor=sum;
        in_site[ini].Seff_nor=100-Sgc-S20;
        in_site[ini].score=r1*in_site[ini].Sspe_nor+r2*in_site[ini].Seff_nor;
        in_site[ini].count=in_site[ini].ot.size();
    }

    int len=in_site[ini].ot.size();
    sort(&(in_site[ini].ot[0]),&(in_site[ini].ot[0])+len,cmp);
    cJSON *otj=Create_array_of_anything(&(in_site[ini].ot[0]),min(20,len));
    sprintf(buffer,"insert into Table_result set result_kind=0, result_gc=%1.4f, result_ntlength=%d, result_Sspe=%.2f, result_Seff=%.2f, result_count=%d, result_offtarget='%s', sgrna_ID=%s; ",in_site[ini].gc,req_restrict.ntlength,in_site[ini].Sspe_nor,in_site[ini].Seff_nor,in_site[ini].count,NomoreSpace(cJSON_Print(otj)),row.row[6]);
    mos_pthread_mutex_lock(&mutex_mysql_conn);
    int res=mysql_query(my_conn,buffer);
    mos_pthread_mutex_unlock(&mutex_mysql_conn);
    if(res){
        printf("%s\n\n",buffer);
        printf("%s\n",mysql_error(my_conn));
        system("pause");
    }
    in_site[ini].otj=otj;
    in_site[ini].ot.clear();
}

/**
@brief Share data between threads.
@see create_thread_socre, new_thread, localrow
*/
struct thread_share_variables{
    localrow *lr;
    localrow row;
    int ini;
    int type;
    double r1;
}thread_share_variables;

/**
@brief New thread function.
*/
void *new_thread(void *args){
    int ini=thread_share_variables.ini;
    score(thread_share_variables.lr,thread_share_variables.row,thread_share_variables.ini,thread_share_variables.type,thread_share_variables.r1);
    mos_pthread_mutex_unlock(&mutex);
    mos_sem_post(&sem_thread);
    in_site[ini].ntid=0;
    mos_pthread_detach(mos_pthread_self());
    return NULL;
}

/**
@brief Function to create new thread.
*/
void create_thread_socre(localrow *lr,localrow row,int ini,int type,double r1){
    mos_pthread_t ntid;
    mos_pthread_mutex_lock(&mutex);
    thread_share_variables.lr=lr;
    thread_share_variables.row=row;
    thread_share_variables.ini=ini;
    thread_share_variables.type=type;
    thread_share_variables.r1=r1;

    int err=mos_pthread_create(&ntid, NULL, new_thread, NULL);
    if(err){
        printf("Pthread_create error: %s\n",strerror(err));
        mos_pthread_mutex_unlock(&mutex);
        mos_sem_post(&sem_thread);
        return ;
    }

    in_site[ini].ntid=ntid;
}
