#include "Main.h"

ptt ptts[1000000];

int pi=0;
site psb_site[1000000];

int ini=0;
site in_site[1000000];

char str[60000000];
char wai[60000000];

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

void print_str(FILE *file,int s,int t){
    int i;
    for(i=s;i<=t;i++){
        fprintf(file,"%c",str[i]);
    }
}

double score(int ini){
    int i;
    int c=0;
    double sum=0;
    for(i=0;i<pi;i++) if(in_site[ini].index!=psb_site[i].index){
        int count=0;
        double ans=0;
        for(int si=0;si<LEN;si++) if(in_site[ini].nt[si]==psb_site[i].nt[si]){
            ans+=M[si];
            count++;
        }
        if(count<20-NUM_NO) ans=0;
        else{
            c++;
            in_site[ini].ot.push_back(i);
        }
        sum+=ans;
    }
    in_site[ini].score=sum;
    in_site[ini].count=c;
    return sum;
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
    // (str[i]=='A' || str[i]=='T' || str[i]=='C' || str[i]=='G') && str[i+1]=='G' && str[i+2]=='G'
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

int main(int args,char *argv[]){
    char ch;

    FILE *ff;
    FILE *fout4;
    FILE *fr=fopen("C:/Users/ZhaoYi/Desktop/NNY-Database/request.txt","r");
    int ri=0;
    char request_str[10000];
    while(fscanf(fr,"%c",&ch)==1){
        request_str[ri++]=ch;
    }
    request_str[ri++]=0;
    cJSON *request=cJSON_Parse(request_str);
    cJSON *cJSON_temp;
    char req_pam[PAM_LEN+1];
    int req_pam_len;
    strcpy(req_pam,cJSON_GetObjectItem(request,"pam")->valuestring);
    req_pam_len=(int)strlen(req_pam);
    char req_specie[50];
    int ptts_num;
    strcpy(req_specie,cJSON_GetObjectItem(request,"specie")->valuestring);
    if(strcmp(req_specie,"SARS")==0){
        ff=fopen("C:/Users/ZhaoYi/Desktop/NNY-Database/SARS.fasta","r");
        fout4=fopen("C:/Users/ZhaoYi/Desktop/NNY-Database/result.out","w");
        ptts_num=ptt_readin(PTT_SARS,ptts);
    }else if(strcmp(req_specie,"E.coli")==0){
        ff=fopen("C:/Users/ZhaoYi/Desktop/NNY-Database/NC_017626.fna","r");
        fout4=fopen("C:/Users/ZhaoYi/Desktop/NNY-Database/result.out","w");
        ptts_num=ptt_readin(PTT_ECOLI,ptts);
    }
    cJSON_temp=cJSON_GetObjectItem(request,"gene");
    int req_gene_start,req_gene_end;
    if(cJSON_temp){
        char req_gene[20];
        strcpy(req_gene,cJSON_temp->valuestring);
        for(int i=0;i<ptts_num;i++){
            if(strcmp(req_gene,ptts[i].gene)==0){
                req_gene_start=ptts[i].s;
                req_gene_end=ptts[i].t;
                break;
            }
        }
    }else{
        char req_location[20];
        strcpy(req_location,cJSON_GetObjectItem(request,"location")->valuestring);
        sscanf(req_location,"%d..%d",&req_gene_start,&req_gene_end);
    }
    /*
    This part above is for read in JSON-style request.
    The result stored in req_specie, req_pam, req_gene_start, req_gene_end and so on.
    At the same time, it reads in ptt file for specie,
    and also open files.
    */

    printf("specie: %s\n",req_specie);
    printf("(start,end): (%d,%d)\n",req_gene_start,req_gene_end);
    printf("pam: %s\n",req_pam);

    int i=1;        //���еĿ�ʼ
    int j=1;    //��ŵĿ�ʼ

    readLine(ff);   //��������
    while(fscanf(ff,"%c",&ch)==1){
        if(ch==10) continue;
        wai[i]=0;
        str[i++]=ch;
    }
    int len=i;      //���г���

    for(i=LEN;i<len-req_pam_len;i++){       // All possible gRNAs, +direction
        if(check_pam(str+i,req_pam)){
            psb_site[pi].index=i;
            psb_site[pi].strand='+';
            strncpy(psb_site[pi].pam,str+i,req_pam_len);
            strncpy(psb_site[pi].nt,str+i-LEN,LEN);
            pi++;
        }
    }
    char req_pam_rev[PAM_LEN];
    dna_rev(req_pam_rev,req_pam,req_pam_len);
    for(i=0;i<len-LEN-req_pam_len;i++){     // All possible gRNAs, -direction
        if(check_pam(str+i,req_pam_rev)){
            psb_site[pi].index=i+req_pam_len-1;
            psb_site[pi].strand='-';
            for(j=0;j<req_pam_len;j++) psb_site[pi].pam[j]=dna_rev_char((str+i)[req_pam_len-1-j]);
            psb_site[pi].pam[j]=0;
            for(j=0;j<LEN;j++) psb_site[pi].nt[j]=dna_rev_char((str+i+req_pam_len)[LEN-j-1]);
            psb_site[pi].nt[j]=0;
            pi++;
        }
    }

    for(j=0;j<ptts_num;j++){        // ���
        int s=ptts[j].s;
        int t=ptts[j].t;
        for(i=s;i<=t;i++){
            wai[i]=1;
        }
    }

    printf("Number of possible gRNA: %d\n",pi);

    for(i=0;i<pi;i++) if(i%10000==0) printf("."); printf("\n");
    for(i=0;i<pi;i++){             //�ҵ� gRNAs ������
        if(i%10000==0) printf(".");
        fflush(stdout);
        if(psb_site[i].strand=='+'){
            if(psb_site[i].index<req_gene_start || psb_site[i].index+req_pam_len-1>req_gene_end) continue;
            for(j=psb_site[i].index-LEN;j<psb_site[i].index+req_pam_len-1;j++){
                if(wai[j]!=1) break;
            }
            if(j<psb_site[i].index+req_pam_len-1) j=0;
            else j=1;
        }else{
            if(psb_site[i].index-req_pam_len+1<req_gene_start || psb_site[i].index>req_gene_end) continue;
            for(j=psb_site[i].index-req_pam_len+1;j<psb_site[i].index+LEN;j++){
                if(wai[j]!=1) break;
            }
            if(j<psb_site[i].index+LEN) j=0;
            else j=1;
        }
        if(j){
            in_site[ini]=psb_site[i];
            in_site[ini].ot.clear();
            score(ini);
            ini++;
        }
    }

    sort(in_site,in_site+ini,cmp_in_site);  // Sort & Output

    cJSON *root=cJSON_CreateObject();
    vector<cJSON*> list;
    list.clear();

    for(i=0;i<ini;i++){
        cJSON *ans=cJSON_CreateObject();
        char buffer[30];
        sprintf(buffer,"#%d",i+1);
        cJSON_AddStringToObject(ans,"key",buffer);
        sprintf(buffer,"%s%s",in_site[i].nt,in_site[i].pam);
        cJSON_AddStringToObject(ans,"grna",buffer);
        cJSON_AddNumberToObject(ans,"position",in_site[i].index);
        cJSON_AddNumberToObject(ans,"total_score",in_site[i].score);
        vector<cJSON*>sublist;
        sublist.clear();
        for(j=0;j<in_site[i].count;j++){
            cJSON *subans=cJSON_CreateObject();
            int x=in_site[i].ot[j];
            sprintf(buffer,"%s%s",psb_site[x].nt,psb_site[x].pam);
            cJSON_AddStringToObject(subans,"osequence",buffer);
            double score=0;
            int omms=LEN;
            for(int si=0;si<20;si++) if(in_site[i].nt[si]==psb_site[x].nt[si]){
                score+=M[si];
                omms--;
            }
            cJSON_AddNumberToObject(subans,"oscore",score);
            cJSON_AddNumberToObject(subans,"omms",omms);
            cJSON_AddNumberToObject(subans,"ostrand",psb_site[x].strand);
            cJSON_AddNumberToObject(subans,"oposition",psb_site[x].index);
            if(psb_site[x].region) strcpy(buffer,"Intergenic");
            else strcpy(buffer,"exco");
            cJSON_AddStringToObject(subans,"oregion",buffer);
            sublist.push_back(subans);
        }
        cJSON_AddItemToObject(ans,"offtarget",Create_array_of_anything(&sublist[0],sublist.size()));
        list.push_back(ans);
    }

    cJSON_AddItemToObject(root,"result",Create_array_of_anything(&list[0],list.size()));

    fprintf(fout4,"%s",cJSON_Print(root));
    return 0;
}
