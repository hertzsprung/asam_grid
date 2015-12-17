#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include "f2c_daten.h"


#ifdef irix64
 fileopenclose_(char *fname, int *fmode, int *ecode)
#else
 fileopenclose(char *fname, int *fmode, int *ecode)
#endif

{
char *zei;
int  i;
/* Test auf Ende String */
   zei=strchr(fname,32);
   i=zei-fname;
   if(i>0) fname[i]='\0';
/*                  */
   *ecode=0;
   switch(*fmode)
       {
         case NEW:{
                 if(!(fp=fopen(fname,"w"))){
                   printf("Fehler beim Kreieren der Datei %s\n",fname);
                   *ecode=11;
                   return(1);
                  }
                 break;
                 }
         case APPEND:{
                 if(!(fp=fopen(fname,"a+"))){
                   printf("Fehler beim Oeffnen der Datei %s\n",fname);
                   *ecode=12;
                   return(1);
                  }
                 break;
                 }
         case READ:{
                 if(!(fp=fopen(fname,"r"))){
                   printf("Fehler beim Oeffnen der Datei %s\n",fname);
                   *ecode=13;
                   return(1);
                  }
                 break;
                 }
         case CLOSE:
                 fclose(fp);
                 fp=NULL;
                 break;
         default:{
                  printf("Unbekannter Filemode");
                  *ecode=13;
                  return(1);
                 }
       }

}

#ifdef irix64
  writef2c_(dtype,dzahl,daten,ecode)
#else
  writef2c(dtype,dzahl,daten,ecode)
#endif

int  *dtype;
int  *dzahl;
void *daten;
int  *ecode; 
{
int    i;

   *ecode=0;
   switch(*dtype)
       {
         case FFLOAT:{
                 fwrite(daten,sizeof(float),*dzahl,fp);
                 break;
                }
         case FDOUBLE:{
                 fwrite(daten,sizeof(double),*dzahl,fp);
                 break;
                }
         case FINT:{
                 fwrite(daten,sizeof(int),*dzahl,fp);
                 break;
                }
         case FLONG :{
                 fwrite(daten,sizeof(long),*dzahl,fp);
                 break;
                }
         case FCHAR:{
                 fwrite(daten,sizeof(char),*dzahl,fp);
                 break;
                }
        default:{
                 printf("Unbekannter Datentyp");
                 *ecode=2;
                 return(1);
                }
       } /*end switch */

    return(0);

}


#ifdef irix64
  readc2f_(dtype,dzahl,daten,ecode)
#else
  readc2f(dtype,dzahl,daten,ecode)
#endif

int  *dtype;
int  *dzahl;
void *daten;
int  *ecode; 

{
   if(feof(fp)) {
     *ecode=-1;
     return(1);
     }

   *ecode=0;
   switch(*dtype)
       {
         case FFLOAT:{
                 fread(daten,sizeof(float),*dzahl,fp);
                 break;
                }
         case FDOUBLE:{
                 fread(daten,sizeof(double),*dzahl,fp);
                 break;
                }
         case FINT:{
                 fread(daten,sizeof(int),*dzahl,fp);
                 break;
                }
         case FLONG:{
                 fread(daten,sizeof(long),*dzahl,fp);
                 break;
                }
         case FCHAR:{
                 fread(daten,sizeof(char),*dzahl,fp);
                 break;
                }
        default:{
                 printf("Unbekannter Datentyp");
                 *ecode=2;
                 return(1);
                }
      } /* end switch */


   return(0);
}

