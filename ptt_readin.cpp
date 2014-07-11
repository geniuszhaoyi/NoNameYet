#include "Main.h"

int ptt_readin(int type,ptt *ptts){
    char buffer1[1000],buffer2[1000];
    int s,t;
    int ptti=0;
    if(type==PTT_SARS){
        FILE *fp=fopen("C:/Users/ZhaoYi/Desktop/SARS/SARS.ptt","r");
        readLine(fp);   //读入外显子区域，并标记
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
    return 0;
}
