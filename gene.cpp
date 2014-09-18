/** @file
@brief Get gene infomation.
@author Yi Zhao
*/
#include "main.h"

/**
@brief Get gene infomation by gene name.
@param str Output, the gene's location.
@param specie_name Input, the specie's name.
@param gene_name Input, the gene's name.
@return 0 for succeed, -1 for failed.
*/
int get_gene_info(char *str,const char *specie_name,const char *gene_name){
    char buffer[1024];
    sprintf(buffer,"SELECT Chr_Name,gene_Start,gene_end FROM Table_gene as g JOIN Table_chromosome as c ON g.Chr_No=c.Chr_No JOIN Table_Specie as s ON c.Sno=s.Sno WHERE SName='%s' and gene_Name='%s'",specie_name,gene_name);
    int res=mysql_query(my_conn,buffer);
    if(res) return -1;
    MYSQL_RES *result=mysql_store_result(my_conn);
    MYSQL_ROW sql_row;
    if((sql_row=mysql_fetch_row(result))){
        sprintf(str,"%s:%s..%s",sql_row[0],sql_row[1],sql_row[2]);
    }
    return 0;
}
