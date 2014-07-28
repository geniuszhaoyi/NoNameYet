#include "main.h"

MYSQL *db_init(){
    my_conn=mysql_init(NULL);
    if(mysql_real_connect(my_conn,"127.0.0.1","root","zy19930108","db",3306,NULL,0))
        return my_conn;
    return NULL;
}

void db_free(){
    mysql_close(my_conn);
}

int main2(void){
    MYSQL *my_conn;
    MYSQL_RES *result;
    MYSQL_ROW sql_row;

    my_conn=mysql_init(NULL);
    if(mysql_real_connect(my_conn,"127.0.0.1","root","zy19930108","db",3306,NULL,0)){
        int res=mysql_query(my_conn,"SELECT * FROM table_specie;");
        if(res)return 1;
        result=mysql_store_result(my_conn);
        while((sql_row=mysql_fetch_row(result))){
            for(int i=0;i<2;i++){
                printf("%s\t",sql_row[i]);
            }
            printf("\n");
        }
    }

    mysql_free_result(result);
    mysql_close(my_conn);

    return 0;
}

