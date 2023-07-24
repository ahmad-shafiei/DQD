#include <stdio.h>

#define MAXLINE 800
#define ZERO 0
#define DIMENSION 330
#define max(a,b) ((a>b) ? a:b)
#define min(a,b) ((a>b) ? b:a)

int numero_isole;
int numero_esterni;
float dCollegamenti[DIMENSION][DIMENSION];

int main(int argc,char **argv)
{
  
  char line[MAXLINE];
  FILE *fp,*out_file;
  int stg,conto=ZERO,caso;

  fp = fopen(argv[1], "r");
  out_file = fopen(argv[2], "w");
  
  inizializzazione(out_file);
  
  while((stg=fgetline(fp,line,MAXLINE)) != EOF){
    conto++;
    switch(caso=line[0]){

    case 35:            /* toglie le righe che iniziano con # */
      break;
      
    case 32:           /* toglie le righe che cominciano con lo spazio */
      if(stg!=0){
	printf("C'e' un errore di formattazione alla riga %d\n",conto);
	exit(1);
      }
      break;
      
    case 84:              /* lettera T */
      scrivi_tunneling(out_file,line);  
      break;

    case 80:              /* lettera P */
      scrivi_parametri(out_file,line,conto);
      break;

    case 67:              /* lettera C */
      scrivi_condensatori(out_file,line);
      break;
      
    case 86:              /* lettera V */
      scrivi_potenziali(out_file,line);
      break;

    default:
      switch(stg){
      case 0:
	break;
	
      default:
	printf("Formato riga %d non riconosciuto",conto);
	exit(1);
      }

    }
  }
  
  chiusura(out_file);

}

/*******************************************************/

inizializzazione(FILE *out_file)
{
  int i,j;

  for(i=0;i<DIMENSION;i++)
    for(j=0;j<DIMENSION;j++)
      dCollegamenti[i][j] = 0;
    fprintf(out_file," &data\n");  
/*  fprintf(out_file," &data\n dLarghezza = 1d-4,\n dSpintina = 1d-4,\n\n"); */
}

chiusura(FILE *out_file)
{
  int i,j;
  for(i=1;i<(numero_isole+1);i++)
    for(j=1;j<(numero_esterni+numero_isole+1);j++)
      fprintf(out_file," dCollegamenti(%d,%d) = %f,\n",
	      i,j,dCollegamenti[i][j]);
  
  for(i=1;i<(numero_isole+1);i++)
    fprintf(out_file," iNn(%d) = %d,\n",i,0);
  
  fprintf(out_file," &end\n");

}
  
scrivi_parametri(FILE *out_file,char line[],int conto)
{
  int numero,caso;
  int valore;
  float temperatura,spintina;
  char descrittore[MAXLINE];
  
  sscanf(line,"%*c%d%s%d",&numero,&descrittore,&valore);

  switch(caso=descrittore[0]){
  case 100:
    sscanf(line,"%*c%*d%*s%f",&spintina);
    fprintf(out_file," dSpintina = %f,\n dLarghezza = %f,\n\n"
        ,spintina,spintina);
    break;  
  case 84:
     sscanf(line,"%*c%*d%*s%f",&temperatura);
    fprintf(out_file," dTemperatura = %f,\n\n",temperatura);
    break;
  case 83:
    switch(caso=descrittore[1]){
    case 101:
      fprintf(out_file," iSeme = %d,\n\n",valore);
      break;
    case 116:
      fprintf(out_file," iPassi = %d,\n\n",valore);
      break;
    default:
      printf("Errore nella riga %d,\n\n",conto);
      exit(1);
    }
    break;
  case 73:
    numero_isole = valore;
    fprintf(out_file," iIsole = %d,\n\n",valore);
    break;
  case 69:
    numero_esterni = valore;
    fprintf(out_file," iNesterni = %d,\n\n",valore);
    break;
  case 74:
    fprintf(out_file," iNtunnel = %d,\n\n",valore);
    break;
  case 65:
    fprintf(out_file," iMedie = %d,\n\n",valore);
    break;
  default:
    printf("Errore nella riga %d,\n\n",conto);
    exit(1);
  }

}

/*******************************************************/


scrivi_condensatori(FILE *out_file,char line[])
{
  int destra, sinistra, minimo, massimo;
  float capacita;
    
  sscanf(line,"%*c %*d %d %d %f %f",
	 &destra,&sinistra,&capacita);

  minimo=min(destra,sinistra);
  massimo=max(destra,sinistra);

  if(massimo <= numero_isole)
    dCollegamenti[massimo][minimo] = capacita;
  dCollegamenti[minimo][massimo] = capacita;


}
/*******************************************************/

scrivi_potenziali(FILE *out_file,char line[])
{
  int numero, sopra, sotto;
  float iniziale, finale;
  
  sscanf(line,"%*c %d %d %d %f %f",&numero,&sopra,&sotto,
	 &iniziale,&finale);
  
  fprintf(out_file,
	  " dVoltaggio(1,%d) = %f,\n dVoltaggio(2,%d) = %f,\n\n",
	  numero,iniziale,numero,finale);
}
/*******************************************************/
scrivi_tunneling(FILE *out_file,char line[])

{
  int quale, destra, sinistra, minimo, massimo;
  float capacita,costante;
    
  sscanf(line,"%*c %d %d %d %f %f",&quale,
	 &destra,&sinistra,&capacita,&costante);

  minimo=min(destra,sinistra);
  massimo=max(destra,sinistra);

  if(massimo <= numero_isole)
    dCollegamenti[massimo][minimo] = capacita;
	    
  dCollegamenti[minimo][massimo] = capacita;
  
  fprintf(out_file," iTunnel(%d,1) = %d,\n iTunnel(%d,2) = %d,\n\n",
	  quale,minimo,quale,massimo);
  
  fprintf(out_file," dCostante(%d) = %f,\n\n", quale,costante);

}
/*******************************************************/

/* Read one line from fp, */
/* copying it to line array (but no more than max chars). */
/* Does not place terminating \n in line array. */
/* Returns line length, or 0 for empty line, or EOF for end-of-file. */

int fgetline(FILE *fp, char line[], int max)
{
int nch = 0;
int c;
max = max - 1;                  /* leave room for '\0' */

while((c = getc(fp)) != EOF)
        {
        if(c == '\n')
                break;

        if(nch < max)
                {
                line[nch] = c;
                nch = nch + 1;
                }
        }

if(c == EOF && nch == 0)
        return EOF;
        
line[nch] = '\0';
return nch;
}



