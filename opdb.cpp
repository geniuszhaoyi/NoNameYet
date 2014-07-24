#include "main.h"

#define K 100000000000000
#define OFFSET 97

/*typedef struct dcset{
    char req_specie[50];
    char req_kind[50];
    char pam[PAM_LEN];
    int type;   //score type
}dcset;*/

bool cmp(cJSON *a,cJSON *b){
    double as=cJSON_GetObjectItem(a,"oscore")->valuedouble;
    double bs=cJSON_GetObjectItem(b,"oscore")->valuedouble;
    return as>bs;
}

void generate_filename(char *str,const char *req_specie,const char *req_kind,const char *pam,const int type){
    *str=0;
    strcpy(str,req_specie);
    strcat(str,";");

    if(req_kind) strcat(str,req_kind);
    strcat(str,";");

    strcat(str,pam);
    strcat(str,";");

    char buffer[3];
    sprintf(buffer,"%d",type);
    strcat(str,buffer);

    strcat(str,".dat");
}

cJSON *dc_root;
map<long long int,int> mp;

char filename[196]="Database/DataCache/";
char buffer[DCFILE_LEN+1];

long long int getIndex(cJSON *node){
    int g;
    long long int p;
    strcpy(buffer,cJSON_GetObjectItem(node,"position")->valuestring);
    sscanf(buffer,"%d:%I64d",&g,&p);
    p+=K*g;
    strcpy(buffer,cJSON_GetObjectItem(node,"strand")->valuestring);
    if(*buffer=='-') p*=-1;
    return p;
}

void dc_init(const char *fn){
    int numread;
    strcat(filename,fn);
    FILE *file;
    file=fopen(filename,"r");
    if(file==NULL){
        FILE *fout=fopen(filename,"w");
        fprintf(fout,"{\"json\":[]}");
        fclose(fout);
        file=fopen(filename,"r");
    }

    numread=fread(buffer,sizeof(char),DCFILE_LEN,file);
    buffer[numread]=0;
    dc_root=cJSON_Parse(buffer);

    fclose(file);

    cJSON *list=cJSON_GetObjectItem(dc_root,"json");
    int len=cJSON_GetArraySize(list);
    for(int i=0;i<len;i++){
        cJSON *node=cJSON_GetArrayItem(list,i);
        mp[getIndex(node)]=i+OFFSET;
    }
}

void dc_save(){
    int numwritten;
    FILE *file=fopen(filename,"w");
    char *p=NomoreSpace(cJSON_Print(dc_root));
    numwritten=fwrite(p,sizeof(char),strlen(p),file);
    if(numwritten==DCFILE_LEN){
        for(numwritten=0;numwritten<25;numwritten++) printf("Error!!!\n");
    }
    fclose(file);
    free(p);
}

cJSON *dc_get(int gene,int position,char strand){
    long long int p=K*gene+position;
    if(strand=='-') p*=-1;
    if(mp[p]==0){
        return NULL;
    }
    int x=mp[p]-OFFSET;
    cJSON *t=cJSON_GetArrayItem(cJSON_GetObjectItem(dc_root,"json"),x);
    if(cJSON_GetObjectItem(t,"islegal")->valueint==0) return NULL;
    return t;
}

cJSON *dc_put(int islegal,int ini){
    cJSON *no=cJSON_CreateObject();
    cJSON_AddNumberToObject(no,"islegal",islegal);

    cJSON_AddStringToObject(no,"sequence",in_site[ini].nt);
    sprintf(buffer,"%d:%d",in_site[ini].chromosome,in_site[ini].index);
    cJSON_AddStringToObject(no,"position",buffer);
    sprintf(buffer,"%c",in_site[ini].strand);
    cJSON_AddStringToObject(no,"strand",buffer);

//    cJSON_AddNumberToObject(no,"score",in_site[ini].score);
    cJSON_AddNumberToObject(no,"Sspe",in_site[ini].Sspe_nor);
    cJSON_AddNumberToObject(no,"Seff",in_site[ini].Seff_nor);

    cJSON_AddNumberToObject(no,"count",in_site[ini].count);

    // get all offtarget
    int i,len=in_site[ini].ot.size();
    cJSON *list[len];
    for(i=0;i<len;i++){
        int x=in_site[ini].ot[i];
        list[i]=cJSON_CreateObject();
        sprintf(buffer,"%s%s",psb_site[x].nt,psb_site[x].pam);
        cJSON_AddStringToObject(list[i],"osequence",buffer);
        int omms;
        double score=subscore(i,x,&omms,0);
        cJSON_AddNumberToObject(list[i],"oscore",score);
        cJSON_AddNumberToObject(list[i],"omms",omms);
        sprintf(buffer,"%d:%d",psb_site[x].chromosome,psb_site[x].index);
        cJSON_AddStringToObject(list[i],"oposition",buffer);
        char xs[2];
        xs[0]=psb_site[x].strand;
        xs[1]=0;
        cJSON_AddStringToObject(list[i],"ostrand",xs);
        if(psb_site[x].region) strcpy(buffer,"Intergenic");
        else strcpy(buffer,"exco");
        cJSON_AddStringToObject(list[i],"oregion",buffer);
    }
    // sort
    sort(list,list+len,cmp);
    // get top 20
    cJSON *otj=Create_array_of_anything(list,min(20,len));
    cJSON_AddItemToObject(no,"offtarget",otj);

    list[0]=cJSON_GetObjectItem(dc_root,"json");
    cJSON_AddItemToArray(list[0],no);
    mp[getIndex(no)]=cJSON_GetArraySize(list[0])-1+OFFSET;

    return otj;
}

/*
#include<stdio.h>
intmain(void)
{
    FILE *stream;
    char list[30];
    int i,numread,numwritten; //Openfileintextmode:
    if((stream=fopen("fread.out","w+t"))!=NULL)
    {
        for(i=0;i<25;i++)
            list[i]=(char)('z'-i);¡¡// Write25characterstostream
        numwritten=fwrite(list,sizeof(char),25,stream);
        printf("Wrote%ditems\n",numwritten);
        fclose(stream);
    }
    else
        printf("Problemopeningthefile\n");
    if((stream=fopen("fread.out","r+t"))!=NULL)
    { //Attempttoreadin25characters
        numread=fread(list,sizeof(char),25,stream);
        printf("Numberofitemsread=%d\n",numread);
        printf("Contentsofbuffer=%.25s\n",list);
        fclose(stream);
    }
    else
        printf("Filecouldnotbeopened\n");
}
*/
