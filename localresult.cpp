#include "main.h"

int make_mysqlres_local(localrow **localresult,MYSQL_RES *result_t){
    int count=0;
    localrow **p=localresult;
    mysql_data_seek(result_t,0);
    MYSQL_ROW sql_row;
    while((sql_row=mysql_fetch_row(result_t))){
        *p=(localrow*)malloc(sizeof(localrow));
        for(int i=0;i<8;i++){
            strcpy((*p)->row[i],sql_row[i]);
        }
        p=&((*p)->next);
        count++;
    }
    *p=NULL;
    return count;
}

void free_mysqlres_local(localrow *localresult){
    while(localresult){
        localrow *p=localresult;
        localresult=p->next;
        free(p);
    }
}

int localres_count(localrow *lr){
    int i=0;
    while(lr){
        i++;
        lr=lr->next;
        if(i>1000000) return -1;
    }
    return i;
}
