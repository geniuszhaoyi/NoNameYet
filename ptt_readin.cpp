#include "Main.h"

int ptt_readin(int type,ptt *ptts){
    char buffer1[1000],buffer2[1000];
    int ptti=0;
    if(type==PTT_SARS){
        int s,t;
        FILE *fp=fopen("C:/Users/ZhaoYi/Desktop/NNY-Database/SARS.ptt","r");
        readLine(fp);   //读入外显子区域
        while(fscanf(fp,"%s%s%d%d",buffer1,buffer2,&s,&t)==4){
            fscanf(fp,"%s%s",buffer1,buffer2);
            readLine(fp);
            ptts[ptti].s=s;
            ptts[ptti].t=t;
            strcpy(ptts[ptti].gene,buffer2);
            ptti++;
        }
        return ptti;
    }
    if(type==PTT_ECOLI){
        int s,t;
        char strand;
        int len,pid;
        FILE *fp=fopen("C:/Users/ZhaoYi/Desktop/NNY-Database/NC_017626.ptt","r");
        readLine(fp);
        readLine(fp);
        readLine(fp);
        while(fscanf(fp,"%d..%d\t%c\t%d\t%d\t%s",&s,&t,&strand,&len,&pid,buffer1)==6){
            readLine(fp);
            ptts[ptti].s=s;
            ptts[ptti].t=t;
            ptts[ptti].strand=strand;
            strcpy(ptts[ptti].gene,buffer1);
            ptti++;
        }
        return ptti;
    }
    return 0;
}
