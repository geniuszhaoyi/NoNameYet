#include "main.h"

struct return_struct info_readin(int type,ptt *ptts,char str[][GENE_LEN],char wai[][GENE_LEN]){
    int ptti=0;
    return_struct rs;
    if(type==PTT_SARS){
        int s,t;
        char buffer1[1000],buffer2[1000];
        FILE *fp=fopen("SARS.ptt","r");
        readLine(fp);   //读入外显子区域
        while(fscanf(fp,"%s%s%d%d",buffer1,buffer2,&s,&t)==4){
            fscanf(fp,"%s%s",buffer1,buffer2);
            readLine(fp);
            ptts[ptti].s=s;
            ptts[ptti].t=t;
            ptts[ptti].chromosome=1;
            strcpy(ptts[ptti].gene,buffer2);
            ptti++;
        }
        rs.ptts_num=ptti;
        return rs;
    }
    if(type==PTT_ECOLI){
        int s,t;
        char strand;
        char buffer[1000];
        int len,pid;
        FILE *fp=fopen("NC_017626.ptt","r");
        readLine(fp);
        readLine(fp);
        readLine(fp);
        while(fscanf(fp,"%d..%d\t%c\t%d\t%d\t%s",&s,&t,&strand,&len,&pid,buffer)==6){
            readLine(fp);
            ptts[ptti].s=s;
            ptts[ptti].t=t;
            ptts[ptti].chromosome=1;
            ptts[ptti].strand=strand;
            strcpy(ptts[ptti].gene,buffer);
            ptti++;
        }
        rs.ptts_num=ptti;
        return rs;
    }
    if(type==PTT_SACCHAROMYCETES){
        FILE *flist=fopen("saccharomycetes/list.txt","r");
        int id;
        int id_max=-1;
        char filename1[64],filename2[64],buffer[128];
        while(fscanf(flist,"%d\t%s\t%s",&id,filename1,filename2)==3){
            int s,t;
            char strand;
            int len,pid;

            strcpy(buffer,"saccharomycetes/");
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

            strcpy(buffer,"saccharomycetes/");
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

            if(id_max<id) id_max=id;
        }
        rs.ptts_num=ptti;
        rs.num_chromosome=id_max;
        return rs;
    }
    return rs;
}
