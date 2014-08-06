#include "main.h"

const double M[]={0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583};
const double eM[]={2.718281828,2.718281828,2.680491036,2.718281828,2.718281828,1.831252209,1.979808257,2.718281828,1.842272751,2.511800935,1.741940985,1.635584119,1.472556491,1.160672989,1.30734714,1.187677833,1.469614321,1.216526905,1.370259311,1.517402513};

bool cmp(cJSON *a,cJSON *b){
    double as=cJSON_GetObjectItem(a,"oscore")->valuedouble;
    double bs=cJSON_GetObjectItem(b,"oscore")->valuedouble;
    return as>bs;
}

int check_region(int i){
    int r=in_site[i].region;
    return req_restrict.region[r];
}

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

//MYSQL_ROW row :  0:sgrna_start, 1:sgrna_end, 2:sgrna_strand, 3:sgrna_seq, 4:sgrna_PAM, 5:Chr_Name, 6:sgrna_ID, 7:Chr_No
cJSON *cJSON_otj(int ini,MYSQL_ROW row,double oscore,int omms){
    char buffer[4096];
    cJSON *root=cJSON_CreateObject();
    sprintf(buffer,"%s%s",row[3],row[4]);
    cJSON_AddStringToObject(root,"osequence",buffer);
    cJSON_AddNumberToObject(root,"oscore",oscore);
    cJSON_AddNumberToObject(root,"omms",omms);
    sprintf(buffer,"%s:%s",row[5],row[0]);
    cJSON_AddStringToObject(root,"oposition",buffer);
    cJSON_AddStringToObject(root,"ostrand",row[2]);
    cJSON_AddStringToObject(root,"oregion",region_info[getRegion(atoi(row[6]),atoi(row[7]),atoi(row[0]),atoi(row[1]))]);
    return root;
}

double subscore(int ini,MYSQL_ROW row,int *Nph,int type){
    int nmm=0;
    int d0=0;
    double smm=0;
    int nph=0;
    int i;
    for(i=0;i<LEN;i++){
        if(in_site[ini].nt[i]!=row[3][i]){
            smm+=eM[i];
            d0+=i;
            nmm++;
        }
    }
    if(nmm==0){
        nph++;
        smm=25.0;
        if(type==1) in_site[ini].ot.push_back(cJSON_otj(ini,row,smm,nmm));
    }else{
        if(nmm<=NUM_NO){
            smm=4.0*smm/(double)nmm/(double)nmm/(4.0*d0/19.0/(double)nmm+1);
            if(type==1) in_site[ini].ot.push_back(cJSON_otj(ini,row,smm,nmm));
        }else{
            smm=0.0;
        }
    }
    if(type==0) (*Nph)=nmm;
    if(type==1) (*Nph)+=nph;
    return smm;
}

void score(MYSQL_RES *result,MYSQL_ROW row,int ini,int type,double r1){
    double r2=1.0-r1;
    int Sgc=0,S20=0;
    int Nph=0;
    double sum=0;
    int gc=0;
    int i;
    char buffer[9182];

    sprintf(buffer,"SELECT sgrna_Sspe, sgrna_Seff, sgrna_count, sgrna_offtarget FROM Table_sgRNA WHERE sgrna_ID=%s and sgrna_offtarget IS NOT NULL",row[6]);
    int res=mysql_query(my_conn,buffer);
    if(res){
    }
    MYSQL_RES *rsdc=mysql_store_result(my_conn);
    MYSQL_ROW sql_row;
    if((sql_row=mysql_fetch_row(rsdc))){
        sscanf(sql_row[0],"%lf",&in_site[ini].Sspe_nor);
        sscanf(sql_row[1],"%lf",&in_site[ini].Seff_nor);
        in_site[ini].score=r1*in_site[ini].Sspe_nor+r2*in_site[ini].Seff_nor;
        sscanf(sql_row[2],"%d",&in_site[ini].count);
        (*pini)++;
        in_site[ini].otj=cJSON_Parse(sql_row[3]);
        return ;
    }

    for(i=0;i<LEN;i++) if(in_site[ini].nt[i]=='C' || in_site[ini].nt[i]=='G') gc++;
    if((double)gc/(double)LEN<0.4 || (double)gc/(double)LEN>0.8) Sgc=65;
    else if((double)gc/(double)LEN>0.5 && (double)gc/(double)LEN<0.7) Sgc=0;
    else Sgc=35;
    if(in_site[ini].nt[19]!='G') S20=35;

    mysql_data_seek(result,0);
    while((sql_row=mysql_fetch_row(result))){
        int start=atoi(sql_row[0]);
        if(in_site[ini].index==start) continue;
        double smm=subscore(ini,sql_row,&Nph,1);
        sum+=smm;
    }

    //sum=sigma+S1
    if(type==1 && Nph>3){
        in_site[ini].Sspe_nor=rs.dou[1]=max(100-sum,0.0);
        in_site[ini].Seff_nor=rs.dou[2]=100-Sgc-S20;
        in_site[ini].score=rs.dou[0]=r1*in_site[ini].Sspe_nor+r2*in_site[ini].Seff_nor;
        in_site[ini].count=in_site[ini].ot.size();
        rs.dou[0]=0.0;
    }else if(type==1){
        in_site[ini].Sspe_nor=rs.dou[1]=max(100-sum,0.0);
        in_site[ini].Seff_nor=rs.dou[2]=100-Sgc-S20;
        in_site[ini].score=rs.dou[0]=r1*in_site[ini].Sspe_nor+r2*in_site[ini].Seff_nor;
        in_site[ini].count=in_site[ini].ot.size();
    }else{
        sum=sum-Sgc-S20+7;
        in_site[ini].Sspe_nor=rs.dou[1]=sum;
        in_site[ini].Seff_nor=rs.dou[2]=100-Sgc-S20;
        in_site[ini].score=rs.dou[0]=r1*in_site[ini].Sspe_nor+r2*in_site[ini].Seff_nor;
        in_site[ini].count=in_site[ini].ot.size();
    }

    int len=in_site[ini].ot.size();
    sort(&(in_site[ini].ot[0]),&(in_site[ini].ot[0])+len,cmp);
    cJSON *otj=Create_array_of_anything(&(in_site[ini].ot[0]),min(20,len));
    sprintf(buffer,"update Table_sgRNA set sgrna_Sspe=%.2f, sgrna_Seff=%.2f, sgrna_count=%d, sgrna_offtarget='%s' where sgrna_ID=%s; ",in_site[ini].Sspe_nor,in_site[ini].Seff_nor,in_site[ini].count,NomoreSpace(cJSON_Print(otj)),row[6]);
    res=mysql_query(my_conn,buffer);
    if(res){
        printf("%s\n\n",buffer);
        printf("%s\n",mysql_error(my_conn));
        system("pause");
    }
    in_site[ini].otj=otj;
}

