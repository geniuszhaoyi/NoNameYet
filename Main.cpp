#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
using namespace std;

#include "cJSON/cJSON.h"

#define LEN 20
#define NUM_NO 7

typedef struct site{
    char nt[LEN+1];
    char pam[3+1];
    int index;
    int count;
    double score;
    vector <int> ot;
}site;

int pi=0;
site psb_site[100000];

int ini=0;
site in_site[100000];

char buffer1[1000],buffer2[1000];
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
            //for(int si=0;si<20+3;si++) printf("%c",ss[si]); printf("\n");
            //printf("%s%s %f\n",psb_site[i].nt,psb_site[i].pam,ans);
        }
        sum+=ans;
    }
    //printf("%f\n\n",sum);
    in_site[ini].score=sum;
    in_site[ini].count=c;
    return sum;
}

int main(int args,char *argv[]){
/*
	char my_json_string[]="{\n\"name\": \"Jack (\\\"Bee\\\") Nimble\", \n\"format\": {\"type\":       \"rect\", \n\"width\":      1920, \n\"height\":     1080, \n\"interlace\":  false,\"frame rate\": 24\n}\n}";
    printf("%s",my_json_string);
    cJSON *root = cJSON_Parse(my_json_string);
    cJSON *format = cJSON_GetObjectItem(root,"format");
    int framerate = cJSON_GetObjectItem(format,"frame rate")->valueint;
    cJSON *n=cJSON_GetObjectItem(format,"framerate");
    printf("\n%d,%d\n",format,n);
    printf("%d\n",framerate);
*/
    return 0;
    FILE *ff=fopen("C:/Users/ZhaoYi/Desktop/SARS/SARS.fasta","r");
    FILE *fp=fopen("C:/Users/ZhaoYi/Desktop/SARS/SARS.ptt","r");
    FILE *fout1=fopen("C:/Users/ZhaoYi/Desktop/SARS/SARS.1.out","w");
    FILE *fout2=fopen("C:/Users/ZhaoYi/Desktop/SARS/SARS.2.out","w");
    FILE *fout3=fopen("C:/Users/ZhaoYi/Desktop/SARS/SARS.3.out","w");
    char ch;
    int i=1;        //序列的开始
    int j=1,k=1;    //编号的开始
    int s,t;

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

//    printf("%s%s %d %d\n",psb_site[0].nt,psb_site[0].pam,psb_site[0].index-LEN,psb_site[0].index+2);
//    printf("%s%s %d %d\n",psb_site[pi-1].nt,psb_site[pi-1].pam,psb_site[pi-1].index-LEN,psb_site[pi-1].index+2);

/*    {
        int i,j;
        int sum=0,n=0;
        for(i=0;i<pi;i++){
            int count=0;
            for(j=0;j<pi;j++){
                int k,c=0;
                for(k=0;k<20;k++){
                    if(psb_site[i].nt[k]==psb_site[j].nt[k]) c++;
                }
                if(c>=20-NUM_NO) count++;
                if(c>=20-NUM_NO && i!=j){
                    printf("%s%s ",psb_site[i].nt,psb_site[i].pam);
                    printf("%s%s\n",psb_site[j].nt,psb_site[i].pam);
                }
            }
            sum+=count;
            n++;
        }
        printf("%d: %d / %d\n",pi,sum,n);
        printf("%f\n",(double)sum/(double)n);
    } */

    readLine(fp);   //读入外显子区域，并标记
    while(fscanf(fp,"%s%s%d%d",buffer1,buffer2,&s,&t)==4){
        readLine(fp);
        fprintf(fout1,"#%d ",j++);
        print_str(fout1,s,t);
        fprintf(fout1," %d %d\n",s,t);
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

    for(i=0;i<ini;i++){
        fprintf(fout3,"#%d %s%s %lf %d %d\n",i+1,in_site[i].nt,in_site[i].pam,in_site[i].score,in_site[i].count,in_site[i].index);
        for(j=0;j<in_site[i].count;j++){
            int x=in_site[i].ot[j];
            fprintf(fout3,"%s%s\n",psb_site[x].nt,psb_site[x].pam);
        }
        fprintf(fout3,"\n");
    }

    /*
    fseek(fp,0,SEEK_SET);
    readLine(fp);
    while(fscanf(fp,"%s%s%d%d",buffer1,buffer2,&s,&t)==4){
        for(i=s;i+2<=t;i++){
            if((str[i]=='A' || str[i]=='T' || str[i]=='C' || str[i]=='G') && str[i+1]=='G' && str[i+2]=='G'){
                print_str(fout2,k++,i-LEN,i+2);
                int ii;
                for(ii=i-LEN;ii<=i+2;ii++){
                    if(wai[ii]==0) break;
                }
                if(ii==i+2) fprintf(fout2," exon\n");
                else fprintf(fout2," inter\n");
            }
        }
        for(i=0;i<pi;i++){
            if(score(psb_site[pi].nt,str+i-20)>=0){
            }
        }
    }
    */
    return 0;
}
