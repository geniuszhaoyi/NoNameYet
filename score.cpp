#include "main.h"

const double M[]={0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583};
const double eM[]={2.718281828,2.718281828,2.680491036,2.718281828,2.718281828,1.831252209,1.979808257,2.718281828,1.842272751,2.511800935,1.741940985,1.635584119,1.472556491,1.160672989,1.30734714,1.187677833,1.469614321,1.216526905,1.370259311,1.517402513};

bool cmp(cJSON *a,cJSON *b){
    double as=cJSON_GetObjectItem(a,"oscore")->valuedouble;
    double bs=cJSON_GetObjectItem(b,"oscore")->valuedouble;
    return as>bs;
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
    cJSON_AddStringToObject(root,"oregion","IDontKnow");
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

void score(MYSQL_RES *result,MYSQL_ROW row,int *pini,int type,double r1){
    double r2=1.0-r1;
    int ini=*pini;
    int Sgc=0,S20=0;
    int Nph=0;
    double sum=0;
    int gc=0;
    int i;

    //MYSQL_ROW row :  sgrna_start, sgrna_end, sgrna_strand, sgrna_seq, sgrna_PAM, Chr_Name
    in_site[ini].index=atoi(row[0]);
    in_site[ini].strand=row[2][0];
    strcpy(in_site[ini].nt,row[3]);
    strcpy(in_site[ini].pam,row[4]);
    in_site[ini].ot.clear();

    return_struct rs;
    if(check_rfc(ini)==0){
        rs.dou[0]=-1.0;
        rs.dou[1]=0.0;
        rs.dou[2]=0.0;
        //dc_put(0,ini);
        return ;
    }

    for(i=0;i<LEN;i++) if(in_site[ini].nt[i]=='C' || in_site[ini].nt[i]=='G') gc++;
    if((double)gc/(double)LEN<0.4 || (double)gc/(double)LEN>0.8) Sgc=65;
    else if((double)gc/(double)LEN>0.5 && (double)gc/(double)LEN<0.7) Sgc=0;
    else Sgc=35;
    if(in_site[ini].nt[19]!='G') S20=35;

    MYSQL_ROW sql_row;
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
        (*pini)++;
        rs.dou[0]=0.0;
    }else if(type==1){
        in_site[ini].Sspe_nor=rs.dou[1]=max(100-sum,0.0);
        in_site[ini].Seff_nor=rs.dou[2]=100-Sgc-S20;
        in_site[ini].score=rs.dou[0]=r1*in_site[ini].Sspe_nor+r2*in_site[ini].Seff_nor;
        in_site[ini].count=in_site[ini].ot.size();
        (*pini)++;
    }else{
        sum=sum-Sgc-S20+7;
        in_site[ini].Sspe_nor=rs.dou[1]=sum;
        in_site[ini].Seff_nor=rs.dou[2]=100-Sgc-S20;
        in_site[ini].score=rs.dou[0]=r1*in_site[ini].Sspe_nor+r2*in_site[ini].Seff_nor;
        in_site[ini].count=in_site[ini].ot.size();
        (*pini)++;
    }

    int len=in_site[ini].ot.size();
    sort(&(in_site[ini].ot[0]),&(in_site[ini].ot[0])+len,cmp);
    cJSON *otj=Create_array_of_anything(&(in_site[ini].ot[0]),min(20,len));
    in_site[ini].otj=otj;
}
