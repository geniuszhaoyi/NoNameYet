#include "Main.h"

ptt ptts[10000];

int pi=0;
site psb_site[100000];

int ini=0;
site in_site[100000];

char str[6000000];
char wai[6000000];

const double M[]={0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583};

bool cmp_in_site(site a,site b){
    return a.score>b.score;
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

double score(char *ss,int ini){
    int i;
    int c=0;
    double sum=0;
    for(i=0;i<pi;i++){
        int count=0;
        double ans=0;
        for(int si=0;si<20;si++) if(ss[si]==psb_site[i].nt[si]){
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

int main(int args,char *argv[]){
    char ch;

    FILE *fr=fopen("C:/Users/ZhaoYi/Desktop/SARS/request.txt","r");
    int ri=0;
    char request_str[10000];
    while(fscanf(fr,"%c",&ch)==1){
        request_str[ri++]=ch;
    }
    request_str[ri++]=0;
    cJSON *request=cJSON_Parse(request_str);
    cJSON *cJSON_temp;
    char req_pam[20];
    strcpy(req_pam,cJSON_GetObjectItem(request,"pam")->valuestring);
    char req_specie[50];
    int ptts_num;
    strcpy(req_specie,cJSON_GetObjectItem(request,"specie")->valuestring);
    if(strcmp(req_specie,"SARS")==0){
        ptts_num=ptt_readin(PTT_SARS,ptts);
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
    The result stored in req_specie, req_gene_start, req_gene_end, req_pam and so on.
    */

    printf("specie: %s\n",req_specie);
    printf("(start,end): (%d,%d)\n",req_gene_start,req_gene_end);
    printf("pam: %s\n",req_pam);

    FILE *ff=fopen("C:/Users/ZhaoYi/Desktop/SARS/SARS.fasta","r");
    FILE *fout2=fopen("C:/Users/ZhaoYi/Desktop/SARS/SARS.2.out","w");
    FILE *fout4=fopen("C:/Users/ZhaoYi/Desktop/SARS/SARS.4.out","w");
    int i=1;        //序列的开始
    int j=1,k=1;    //编号的开始

    readLine(ff);   //读入序列
    while(fscanf(ff,"%c",&ch)==1){
        if(ch==10) continue;
        wai[i]=0;
        str[i++]=ch;
    }
    int len=i;      //序列长度

    for(i=LEN+1;i<len;i++){
        if((str[i]=='A' || str[i]=='T' || str[i]=='C' || str[i]=='G') && (str[i+1]=='G' || str[i+1]=='A') && str[i+2]=='G'){
            psb_site[pi].index=i;
            strncpy(psb_site[pi].pam,str+i,3);
            strncpy(psb_site[pi].nt,str+i-LEN,LEN);
            pi++;
        }
    }

    for(j=0;j<ptts_num;j++){
        int s=ptts[j].s;
        int t=ptts[j].t;
        for(i=s;i<=t;i++){
            wai[i]=1;
        }
    }

    for(i=20;i<len;i++){        //正向找到 NGG
        if(wai[i]==1 && wai[i+1]==1 && wai[i+2]==1){
            if((str[i]=='A' || str[i]=='T' || str[i]=='C' || str[i]=='G') && str[i+1]=='G' && str[i+2]=='G'){
                in_site[ini].index=i;
                strncpy(in_site[ini].pam,str+i,3);
                strncpy(in_site[ini].nt,str+i-LEN,LEN);
                in_site[ini].ot.clear();
                score(str+i-LEN,ini);
                ini++;

                fprintf(fout2,"#%d ",k++);
                print_str(fout2,i-LEN,i+2);
                int ii;
                for(ii=i-LEN;ii<=i+2;ii++){
                    if(wai[ii]==0) break;
                }
                fprintf(fout2," %d %d ",i-LEN,i+2);
                if(ii==i+3) fprintf(fout2,"exon\n");
                else fprintf(fout2,"inter\n");
            }
        }
    }

    sort(in_site,in_site+ini,cmp_in_site);

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
