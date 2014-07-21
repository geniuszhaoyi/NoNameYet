#include "main.h"

struct return_struct info_readin(int type,ptt *ptts,char str[][GENE_LEN],char wai[][GENE_LEN],const char *ctrl){
    int ptti=0;
    return_struct rs;

    if(type==PTT_SARS){
        char path[196];
        char Path[196]="Database/SARS/";
        sprintf(path,"%slist.txt",Path);
        FILE *flist=fopen(path,"r");

        int id;
        int id_max=-1;
        char filename1[64],filename2[64];
        while(fscanf(flist,"%d\t%s\t%s",&id,filename1,filename2)==3){
            int i=1;
            int s,t;
            char ch;
            char buffer1[1000],buffer2[1000];

            sprintf(path,"%s%s",Path,filename1);
            FILE *ff=fopen(path,"r");
            readLine(ff);
            while(fscanf(ff,"%c",&ch)==1){
                if(ch==10) continue;
                wai[id][i]=0;
                str[id][i++]=ch;
            }
            rs.len[id]=i;
            fclose(ff);

            sprintf(path,"%s%s",Path,filename2);
            FILE *fp=fopen(path,"r");
            readLine(fp);
            while(fscanf(fp,"%s%s%d%d",buffer1,buffer2,&s,&t)==4){
                readLine(fp);
                ptts[ptti].s=s;
                ptts[ptti].t=t;
                ptts[ptti].chromosome=id;
                ptts[ptti].strand='.';
                strcpy(ptts[ptti].gene,buffer2);
                ptti++;
                for(i=s;i<=t;i++){
                    wai[id][i]=1;
                }
            }
            fclose(fp);

            if(id_max<id) id_max=id;
        }

        rs.ptts_num=ptti;
        rs.num_chromosome=id_max;
        return rs;
    }

    if(type==PTT_ECOLI){
        char path[196];
        char Path[196]="Database/E.coli/";
        strcat(Path,ctrl);
        strcat(Path,"/");
        sprintf(path,"%slist.txt",Path);
        FILE *flist=fopen(path,"r");

        int id;
        int id_max=-1;
        char filename1[64],filename2[64];
        while(fscanf(flist,"%d\t%s\t%s",&id,filename1,filename2)==3){
            int i=1;
            int s,t;
            char strand,ch;
            int len,pid;
            char buffer[1000];

            sprintf(path,"%s%s",Path,filename1);
            FILE *ff=fopen(path,"r");
            readLine(ff);
            while(fscanf(ff,"%c",&ch)==1){
                if(ch==10) continue;
                wai[id][i]=0;
                str[id][i++]=ch;
            }
            rs.len[id]=i;
            fclose(ff);

            sprintf(path,"%s%s",Path,filename2);
            FILE *fp=fopen(path,"r");
            readLine(fp);
            readLine(fp);
            readLine(fp);
            while(fscanf(fp,"%d..%d\t%c\t%d\t%d\t%s",&s,&t,&strand,&len,&pid,buffer)==6){
                readLine(fp);
                ptts[ptti].s=s;
                ptts[ptti].t=t;
                ptts[ptti].chromosome=id;
                ptts[ptti].strand=strand;
                strcpy(ptts[ptti].gene,buffer);
                ptti++;
                for(i=s;i<=t;i++){
                    wai[id][i]=1;
                }
            }
            fclose(fp);

            if(id_max<id) id_max=id;
        }
        rs.ptts_num=ptti;
        rs.num_chromosome=id_max;
        return rs;
    }

    if(type==PTT_SACCHAROMYCETES){
        FILE *flist=fopen("Database/saccharomycetes/list.txt","r");
        int id;
        int id_max=-1;
        char filename1[64],filename2[64],buffer[128];
        while(fscanf(flist,"%d\t%s\t%s",&id,filename1,filename2)==3){
            int s,t;
            char strand;
            int len,pid;

            strcpy(buffer,"Database/saccharomycetes/");
            strcat(buffer,filename1);
            FILE *ff=fopen(buffer,"r");
            char ch;
            int i=1;
            readLine(ff);
            while(fscanf(ff,"%c",&ch)==1){
                if(ch==10) continue;
                wai[id][i]=0;
                str[id][i++]=ch;
            }
            rs.len[id]=i;
            fclose(ff);

            strcpy(buffer,"Database/saccharomycetes/");
            strcat(buffer,filename2);
            FILE *fp=fopen(buffer,"r");
            readLine(fp);
            readLine(fp);
            readLine(fp);
            while(fscanf(fp,"%d..%d\t%c\t%d\t%d\t%s",&s,&t,&strand,&len,&pid,buffer)==6){
                readLine(fp);
                ptts[ptti].s=s;
                ptts[ptti].t=t;
                ptts[ptti].chromosome=id;
                ptts[ptti].strand=strand;
                strcpy(ptts[ptti].gene,buffer);
                ptti++;
                for(i=s;i<=t;i++){
                    wai[id][i]=1;
                }
            }
            fclose(fp);

            if(id_max<id) id_max=id;
        }
        rs.ptts_num=ptti;
        rs.num_chromosome=id_max;
        return rs;
    }
    return rs;
}
