#include "main.h"

//req_GobalParameter
restrict req_restrict;

ptt ptts[1000000];

int pi;
site psb_site[1000000];

int ini;
site in_site[1000000];

char str[NUM_CHROMOSOME][GENE_LEN];
char wai[NUM_CHROMOSOME][GENE_LEN];

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

char argv_default[]="{\"specie\":\"E.coli\",\"kind\":\"E.coli K12-MG1655\",\"location\":\"1:1..500\",\"pam\":\"NGG\",\"rfc\":\"100010\"}";

int main(int args,char *argv[]){
    const int smallOutputNumber=-1;
    int i,j;
    cJSON *root;
    int num_chromosome;
    int req_type=1;
    char buffer[30];
    char req_gene[20]={0};
    cJSON *msg;
    char req_kind[50]="E.coli K12-DH10B";

    pi=0;
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
    j=1;

    cJSON *cJSON_temp;
    cJSON_temp=cJSON_GetObjectItem(request,"type");
    if(cJSON_temp) req_type=cJSON_temp->valueint;
    char req_pam[PAM_LEN+1];
    int req_pam_len;
    strcpy(req_pam,cJSON_GetObjectItem(request,"pam")->valuestring);
    req_pam_len=(int)strlen(req_pam);
    char req_specie[30];
    struct return_struct rs;
    int ptts_num;
    int len[NUM_CHROMOSOME];
    strcpy(req_specie,cJSON_GetObjectItem(request,"specie")->valuestring);
    if(strcmp(req_specie,"SARS")==0){
        rs=info_readin(PTT_SARS,ptts,str,wai,NULL);
    }else if(strcmp(req_specie,"E.coli")==0){
        cJSON_temp=cJSON_GetObjectItem(request,"kind");
        if(cJSON_temp) strcpy(req_kind,cJSON_temp->valuestring);
        rs=info_readin(PTT_ECOLI,ptts,str,wai,req_kind);
    }else if(strcmp(req_specie,"Saccharomycetes")==0){
        rs=info_readin(PTT_SACCHAROMYCETES,ptts,str,wai,NULL);
    }else{
        onError("no specie");
        return 0;
    }
    ptts_num=rs.ptts_num;
    num_chromosome=rs.num_chromosome;
    for(i=1;i<=num_chromosome;i++) len[i]=rs.len[i];

    double req_r1=0.65;
    cJSON_temp=cJSON_GetObjectItem(request,"gene");
    if(cJSON_temp){
        req_r1=cJSON_GetObjectItem(request,"r1")->valuedouble;
    }

    cJSON_temp=cJSON_GetObjectItem(request,"gene");
    int req_id,req_gene_start,req_gene_end;
    if(cJSON_temp){
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

    generate_filename(buffer,req_specie,req_kind,req_pam,req_type);
    dc_init(buffer);
    /*
    This part above is for read in JSON-style request.
    The result stored in req_specie, req_pam, req_gene_start, req_gene_end, rfc and so on.
    At the same time, it reads in ptt file for specie,
    and also open files.
    */

    for(int id=1;id<=num_chromosome;id++){
        for(i=LEN;i<len[id]-req_pam_len;i++){       // All possible gRNAs, +direction
            if(check_pam(str[id]+i,req_pam)){
                psb_site[pi].index=i;
                psb_site[pi].strand='+';
                psb_site[pi].chromosome=id;
                for(j=0;j<req_pam_len;j++) psb_site[pi].pam[j]=(str[id]+i)[j];
                psb_site[pi].pam[j]=0;
                for(j=0;j<LEN;j++) psb_site[pi].nt[j]=(str[id]+i-LEN)[j];
                psb_site[pi].nt[j]=0;
                pi++;
            }
        }
        char req_pam_rev[PAM_LEN];
        int last=-1;
        dna_rev(req_pam_rev,req_pam,req_pam_len);
        for(i=0;i<len[id]-LEN-req_pam_len;i++){     // All possible gRNAs, -direction
            if(check_pam(str[id]+i,req_pam_rev)){
                psb_site[pi].index=i+req_pam_len-1;
                psb_site[pi].strand='-';
                psb_site[pi].chromosome=id;
                if(req_pam_len!=last){
                    last=req_pam_len;
                }
                for(j=0;j<req_pam_len;j++) psb_site[pi].pam[j]=dna_rev_char((str[id]+i)[req_pam_len-1-j]);
                psb_site[pi].pam[j]=0;
                for(j=0;j<LEN;j++) psb_site[pi].nt[j]=dna_rev_char((str[id]+i+req_pam_len)[LEN-j-1]);
                psb_site[pi].nt[j]=0;
                pi++;
            }
        }
    }

    for(i=0;i<pi;i++){
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
            score(i,&ini,req_type,req_r1);
        }
    }

    dc_save();

    sort(in_site,in_site+ini,cmp_in_site);  // Sort & Output

    root=cJSON_CreateObject();
    cJSON_AddNumberToObject(root,"status",0);

    msg=cJSON_CreateObject();
    cJSON_AddStringToObject(msg,"specie",req_specie);
    cJSON_AddStringToObject(msg,"kind",req_kind);
    cJSON_AddStringToObject(msg,"gene",req_gene);
    sprintf(buffer,"%d:%d..%d",req_id,req_gene_start,req_gene_end);
    cJSON_AddStringToObject(msg,"location",buffer);
    cJSON_AddItemToObject(root,"message",msg);

    vector<cJSON*> list;
    list.clear();
    for(i=0;i<ini && i!=smallOutputNumber;i++){
        cJSON *ans=cJSON_CreateObject();
        sprintf(buffer,"#%d",i+1);
        cJSON_AddStringToObject(ans,"key",buffer);
        sprintf(buffer,"%s%s",in_site[i].nt,in_site[i].pam);
        cJSON_AddStringToObject(ans,"grna",buffer);
        sprintf(buffer,"%d:%d",in_site[i].chromosome,in_site[i].index);
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

    printf("%s\n",NomoreSpace(argv[0]=cJSON_Print(root)));

    free(argv[0]);
    return 0;
}
