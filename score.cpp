#include "main.h"

const double M[]={0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583};
const double eM[]={1,1,1.014098459,1,1,1.484384191,1.373002572,1,1.475504551,1.082204322,1.560490196,1.661963941,1.845960983,2.341987669,2.079234922,2.288736686,1.8496566,2.23446092,1.983771836,1.791404591};

int check_rfc(int i){
    char str[LEN+PAM_LEN+3];
    strcpy(str,psb_site[i].nt);
    strcat(str,psb_site[i].pam);
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

double subscore(int ini,int j,int *Nph,int type){
    int nmm=0;
    int d0=0;
    double smm=0;
    int nph=0;
    int i;
    for(i=0;i<LEN;i++){
        if(in_site[ini].nt[i]!=psb_site[j].nt[i]){
            smm+=eM[i];
            d0+=19-i;
            nmm++;
        }
    }
    if(nmm==0){
        nph++;
        if(type==1) in_site[ini].ot.push_back(j);
        smm=25.0;
    }else{
        if(nmm<=NUM_NO){
            smm=2.0*smm/(double)nmm/(double)nmm/(4.0*d0/19.0/(double)nmm+1);
            if(type==1) in_site[ini].ot.push_back(j);
        }else{
            smm=0.0;
        }
    }
    if(type==0) (*Nph)=nmm;
    if(type==1) (*Nph)+=nph;
    return smm;
}

return_struct score(int ii,int *pini,int type,double r1){
    double r2=1.0-r1;
    int ini=*pini;
    int Sgc=0,S20=0;
    int Nph=0;
    double sum=0;
    int gc=0;
    int i;
    return_struct rs;

    in_site[ini]=psb_site[ii];
    in_site[ini].ot.clear();

    if(check_rfc(ini)==0){
        rs.dou[0]=-1.0;
        rs.dou[1]=0.0;
        rs.dou[2]=0.0;
        dc_put(0,ini);
        return rs;
    }

    cJSON *cache=dc_get(in_site[ini].chromosome,in_site[ini].index,in_site[ini].strand);
    if(cache!=NULL){
        //in_site[ini].score=cJSON_GetObjectItem(cache,"score")->valuedouble;
        in_site[ini].Sspe_nor=rs.dou[1]=cJSON_GetObjectItem(cache,"Sspe")->valuedouble;
        in_site[ini].Seff_nor=rs.dou[2]=cJSON_GetObjectItem(cache,"Seff")->valuedouble;
        in_site[ini].score=rs.dou[0]=r1*in_site[ini].Sspe_nor+r2*in_site[ini].Seff_nor;
        in_site[ini].count=cJSON_GetObjectItem(cache,"count")->valuedouble;
        in_site[ini].otj=cJSON_GetObjectItem(cache,"offtarget");
        (*pini)++;
        return rs;
    }

    for(i=0;i<LEN;i++) if(in_site[ini].nt[i]=='C' || in_site[ini].nt[i]=='G') gc++;
    if((double)gc/(double)LEN<0.4 || (double)gc/(double)LEN>0.8) Sgc=65;
    else if((double)gc/(double)LEN>0.5 && (double)gc/(double)LEN<0.7) Sgc=0;
    else Sgc=35;
    if(in_site[ini].nt[19]!='G') S20=35;

    for(int j=0;j<pi;j++) if(in_site[ini].index!=psb_site[j].index){
        double smm=subscore(ini,j,&Nph,1);
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

    in_site[ini].otj=dc_put(1,ini);

    return rs;
}
