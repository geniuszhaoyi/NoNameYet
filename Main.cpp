#define EXPORT_HELLO_DLL
#include "Main.h"

restrict req_restrict;

ptt ptts[1000000];

int pi;
site psb_site[1000000];

int ini;
site in_site[1000000];

char str[NUM_CHROMOSOME][GENE_LEN];
char wai[NUM_CHROMOSOME][GENE_LEN];

const double M[]={0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583};

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

int check_pam(const char *str,const char *pam){
    for(;*pam;pam++,str++){
        if(*pam=='N' || *pam=='n') continue;
        if(*str!=*pam) return 0;
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
            smm+=M[i];
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
            smm=10.0*smm/(double)nmm/(double)nmm/(4.0*d0/19.0/(double)nmm+1);
            if(type==1) in_site[ini].ot.push_back(j);
        }else{
            smm=0.0;
        }
    }
    if(type==0) (*Nph)=nmm;
    if(type==1) (*Nph)+=nph;
    return smm;
}

double score(int ii,int *pini){
    int ini=*pini;
    int Sgc=0,S20=0;
    int Nph=0;
    double sum=0;
    int gc=0;
    int i;

    if(check_rfc(ini)==0) return -1.0;

    in_site[ini]=psb_site[ii];
    in_site[ini].ot.clear();

    if(in_site[ini].index==2639)
        i=100;

    for(i=0;i<LEN;i++) if(in_site[ini].nt[i]=='C' || in_site[ini].nt[i]=='G') gc++;
    if((double)gc/(double)LEN<0.4 || (double)gc/(double)LEN>0.8) Sgc=5;
    if(in_site[ini].nt[19]!='G') S20=2;

    for(int j=0;j<pi;j++) if(in_site[ini].index!=psb_site[j].index){
        double smm=subscore(ini,j,&Nph,1);
        sum+=smm;
    }
    if(Nph>3){
        in_site[ini].score=0.0;
        in_site[ini].count=in_site[ini].ot.size();
        (*pini)++;
        return 0.0;
    }else{
        sum=100-sum-Sgc-S20;
        in_site[ini].score=sum;
        in_site[ini].count=in_site[ini].ot.size();
        (*pini)++;
        return sum;
    }
}

HELLO_API char *test(char *argv,int smallOutputNumber){
    int i,j;

    pi=0;
    ini=0;
    req_restrict.rfc10=0;
    req_restrict.rfc12=0;
    req_restrict.rfc12a=0;
    req_restrict.rfc21=0;
    req_restrict.rfc23=0;
    req_restrict.rfc25=0;
    for(i=0;i<NUM_CHROMOSOME;i++) for(j=0;j<GENE_LEN;j++){
        wai[i][j]=0;
        str[i][j]=0;
    }

    i=1;        //序列的开始
    j=1;        //编号的开num_chromosome始

    cJSON *request=cJSON_Parse(argv);
    cJSON *cJSON_temp;
    char req_pam[PAM_LEN+1];
    int req_pam_len;
    strcpy(req_pam,cJSON_GetObjectItem(request,"pam")->valuestring);
    req_pam_len=(int)strlen(req_pam);
    char req_specie[50];

    struct return_struct rs;
    int ptts_num;
    int len[NUM_CHROMOSOME];
    strcpy(req_specie,cJSON_GetObjectItem(request,"specie")->valuestring);
    if(strcmp(req_specie,"SARS")==0){
        rs=info_readin(PTT_SARS,ptts,str,wai,NULL);
    }else if(strcmp(req_specie,"E.coli")==0){
        cJSON_temp=cJSON_GetObjectItem(request,"kind");
        if(cJSON_temp) rs=info_readin(PTT_ECOLI,ptts,str,wai,cJSON_temp->valuestring);
        else rs=info_readin(PTT_ECOLI,ptts,str,wai,NULL);
    }else if(strcmp(req_specie,"Saccharomycetes")==0){
        rs=info_readin(PTT_SACCHAROMYCETES,ptts,str,wai,NULL);
    }
    ptts_num=rs.ptts_num;
    int num_chromosome=rs.num_chromosome;
    for(i=1;i<=num_chromosome;i++) len[i]=rs.len[i];

    cJSON_temp=cJSON_GetObjectItem(request,"gene");
    int req_id,req_gene_start,req_gene_end;
    if(cJSON_temp){
        char req_gene[20];
        strcpy(req_gene,cJSON_temp->valuestring);
        for(int i=0;i<ptts_num;i++){
            if(strcmp(req_gene,ptts[i].gene)==0){
                req_id=ptts[i].chromosome;
                req_gene_start=ptts[i].s;
                req_gene_end=ptts[i].t;
                break;
            }
        }
    }else{
        char req_location[20];
        strcpy(req_location,cJSON_GetObjectItem(request,"location")->valuestring);
        sscanf(req_location,"%d:%d..%d",&req_id,&req_gene_start,&req_gene_end);
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

    printf("specie: %s\n",req_specie);
    printf("(start,end): (%d,%d)\n",req_gene_start,req_gene_end);
    printf("pam: %s\n",req_pam);

    for(int id=1;id<=num_chromosome;id++){
        for(i=LEN;i<len[id]-req_pam_len;i++){       // All possible gRNAs, +direction
            if(check_pam(str[id]+i,req_pam)){
                psb_site[pi].index=i;
                psb_site[pi].strand='+';
                psb_site[pi].chromosome=id;
                strncpy(psb_site[pi].pam,str[id]+i,req_pam_len);
                strncpy(psb_site[pi].nt,str[id]+i-LEN,LEN);
                pi++;
            }
        }
        char req_pam_rev[PAM_LEN];
        dna_rev(req_pam_rev,req_pam,req_pam_len);
        for(i=0;i<len[id]-LEN-req_pam_len;i++){     // All possible gRNAs, -direction
            if(check_pam(str[id]+i,req_pam_rev)){
                psb_site[pi].index=i+req_pam_len-1;
                psb_site[pi].strand='-';
                psb_site[pi].chromosome=id;
                for(j=0;j<req_pam_len;j++) psb_site[pi].pam[j]=dna_rev_char((str[id]+i)[req_pam_len-1-j]);
                psb_site[pi].pam[j]=0;
                for(j=0;j<LEN;j++) psb_site[pi].nt[j]=dna_rev_char((str[id]+i+req_pam_len)[LEN-j-1]);
                psb_site[pi].nt[j]=0;
                pi++;
            }
        }
    }

    printf("Number of possible gRNA: %d\n",pi);

    for(i=0;i<pi;i++) if(i%20000==0) printf("."); printf("\n");
    for(i=0;i<pi;i++){             //找到 gRNAs 并评分
        if(i%20000==0) printf(".");
        fflush(stdout);
        if(psb_site[i].chromosome!=req_id) continue;
        if(psb_site[i].strand=='+'){
            if(psb_site[i].index<req_gene_start || psb_site[i].index+req_pam_len-1>req_gene_end) continue;
            for(j=psb_site[i].index-LEN;j<psb_site[i].index+req_pam_len-1;j++){
                if(wai[req_id][j]!=1) break;
            }
            if(j<psb_site[i].index+req_pam_len-1) j=0;
            else j=1;
        }else{
            if(psb_site[i].index-req_pam_len+1<req_gene_start || psb_site[i].index>req_gene_end) continue;
            for(j=psb_site[i].index-req_pam_len+1;j<psb_site[i].index+LEN;j++){
                if(wai[req_id][j]!=1) break;
            }
            if(j<psb_site[i].index+LEN) j=0;
            else j=1;
        }
        if(j){
            score(i,&ini);
        }
    }

    sort(in_site,in_site+ini,cmp_in_site);  // Sort & Output

    cJSON *root=cJSON_CreateObject();
    vector<cJSON*> list;
    list.clear();

    for(i=0;i<ini && i!=smallOutputNumber;i++){
        cJSON *ans=cJSON_CreateObject();
        char buffer[30];
        sprintf(buffer,"#%d",i+1);
        cJSON_AddStringToObject(ans,"key",buffer);
        sprintf(buffer,"%s%s",in_site[i].nt,in_site[i].pam);
        cJSON_AddStringToObject(ans,"grna",buffer);
        sprintf(buffer,"%d:%d",in_site[i].chromosome,in_site[i].index);
        cJSON_AddStringToObject(ans,"position",buffer);
        cJSON_AddNumberToObject(ans,"total_score",in_site[i].score);
        cJSON_AddNumberToObject(ans,"count",in_site[i].count);
        vector<cJSON*>sublist;
        sublist.clear();
        for(j=0;j<in_site[i].count && j!=smallOutputNumber;j++){
            cJSON *subans=cJSON_CreateObject();
            int x=in_site[i].ot[j];
            sprintf(buffer,"%s%s",psb_site[x].nt,psb_site[x].pam);
            cJSON_AddStringToObject(subans,"osequence",buffer);
            int omms;
            double score=subscore(i,x,&omms,0);
            cJSON_AddNumberToObject(subans,"oscore",score);
            cJSON_AddNumberToObject(subans,"omms",omms);
            cJSON_AddNumberToObject(subans,"ostrand",psb_site[x].strand);
            sprintf(buffer,"%d:%d",psb_site[x].chromosome,psb_site[x].index);
            cJSON_AddStringToObject(subans,"oposition",buffer);
            if(psb_site[x].region) strcpy(buffer,"Intergenic");
            else strcpy(buffer,"exco");
            cJSON_AddStringToObject(subans,"oregion",buffer);
            sublist.push_back(subans);
        }
        cJSON_AddItemToObject(ans,"offtarget",Create_array_of_anything(&sublist[0],sublist.size()));
        list.push_back(ans);
    }

    cJSON_AddItemToObject(root,"result",Create_array_of_anything(&list[0],list.size()));

    return cJSON_Print(root);
}
