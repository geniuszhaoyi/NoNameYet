#include "main.h"

const double M[]={0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583};
const double eM[]={2.718281828,2.718281828,2.680491036,2.718281828,2.718281828,1.831252209,1.979808257,2.718281828,1.842272751,2.511800935,1.741940985,1.635584119,1.472556491,1.160672989,1.30734714,1.187677833,1.469614321,1.216526905,1.370259311,1.517402513};

bool cmp(cJSON *a,cJSON *b){
    double as=cJSON_GetObjectItem(a,"oscore")->valuedouble;
    double bs=cJSON_GetObjectItem(b,"oscore")->valuedouble;
    return as>bs;
}

//MYSQL_ROW row :  0:sgrna_start, 1:sgrna_end, 2:sgrna_strand, 3:sgrna_seq, 4:sgrna_PAM, 5:Chr_Name, 6:sgrna_ID, 7:Chr_No
cJSON *cJSON_otj(int ini,localrow *lr,double oscore,int omms){
    char buffer[4096];
    cJSON *root=cJSON_CreateObject();
    sprintf(buffer,"%s%s",lr->row[3],lr->row[4]);
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
        if(type==1) in_site[ini].ot.push_back(cJSON_otj(ini,lr,smm,nmm));
    }else{
        if(nmm<=NUM_NO){
            smm=4.0*smm/(double)nmm/(double)nmm/(4.0*d0/((double)req_restrict.ntlength-1.0)/(double)nmm+1);
            if(type==1) in_site[ini].ot.push_back(cJSON_otj(ini,lr,smm,nmm));
        }else{
            smm=0.0;
        }
    }
    if(type==0) (*Nph)=nmm;
    if(type==1) (*Nph)+=nph;
    return smm;
}

void score(localrow *lr,localrow row,int ini,int type,double r1){
    double r2=1.0-r1;
    int Sgc=0,S20=0;
    int Nph=0;
    double sum=0;
    int gc=0;
    int i;
    char buffer[9182];

    mos_pthread_mutex_unlock(&mutex);
    for(i=0;i<LEN;i++) if(in_site[ini].nt[i]=='C' || in_site[ini].nt[i]=='G') gc++;
    if((double)gc/(double)LEN<0.4 || (double)gc/(double)LEN>0.8) Sgc=65;
    else if((double)gc/(double)LEN>0.5 && (double)gc/(double)LEN<0.7) Sgc=0;
    else Sgc=35;
    if(in_site[ini].nt[LEN-req_restrict.ntlength]!='G') S20=35;

    while(lr){
    //printf("12.1\n");
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
    //sprintf(buffer,"update Table_sgRNA set sgrna_Sspe=%.2f, sgrna_Seff=%.2f, sgrna_count=%d, sgrna_offtarget='%s' where sgrna_ID=%s; ",in_site[ini].Sspe_nor,in_site[ini].Seff_nor,in_site[ini].count,NomoreSpace(cJSON_Print(otj)),row.row[6]);
    sprintf(buffer,"insert into Table_result set result_kind=0, result_ntlength=%d, result_Sspe=%.2f, result_Seff=%.2f, result_count=%d, result_offtarget='%s', sgrna_ID=%s; ",req_restrict.ntlength,in_site[ini].Sspe_nor,in_site[ini].Seff_nor,in_site[ini].count,NomoreSpace(cJSON_Print(otj)),row.row[6]);
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

struct thread_share_variables{
    localrow *lr;
    localrow row;
    int ini;
    int type;
    double r1;
}thread_share_variables;

void *new_thread(void *args){
    int ini=thread_share_variables.ini;
    score(thread_share_variables.lr,thread_share_variables.row,thread_share_variables.ini,thread_share_variables.type,thread_share_variables.r1);
    mos_pthread_mutex_unlock(&mutex);
    mos_sem_post(&sem_thread);
    in_site[ini].ntid=0;
    mos_pthread_detach(mos_pthread_self());
    return NULL;
}

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
