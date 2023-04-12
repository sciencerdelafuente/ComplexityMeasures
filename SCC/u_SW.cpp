/**************************************************************************************/  
//  Computes the Sequence Compositional Complexity (SCC_SW) of a genome in fasta format
//  (https://link.aps.org/doi/10.1103/PhysRevLett.80.1344)
//  Related works:  -> https://doi.org/10.3390/biology12020322 
//                  -> https://doi.org/10.1038/s41598-020-76014-4
/**************************************************************************************/  
//  Code written by Rebeca de la Fuente :  2022                                 
/**************************************************************************************/  
//** Some Author Notes:
//** SCC_SW con significaciones desde 0.99 a 0.90                     **\\
// S[] es un vector compuesto por carácteres ACGT de forma ordenada, y no contiene los carácteres 'N'
//Cortes[] es un vector que guarda los cortes realizados. Si hay cadenas de 'N' de longitud mayor a 10,
//S[] registra un corte que separa el genoma por éstas cadenas largas de 'N' presentes en el archivo.

//** Este programa saca por pantalla:
//  >> Tamaño genoma  ---  carácteres totales del archivo
//  >> Cortes iniciales correspondientes a las cadenas de nucleótidos desconocidos,'N'
//  >> Entropía máxima 
//  >> JS inicial (con respecto a los cortes iniciales), que será restado a cada resultado final.
//  >> Nivel de significación impuesta - SCC - SCC*sig - número de cortes          **\\
//  >> Nivel de sig para la cual la SCC se maximiza  -  SCCmax

#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<cstdlib>
#include<string.h>
#include<fstream>
#include<time.h>
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30


using namespace std;

float sl_m_d(int dim,float d, int n);
float gammp_n_f(float x, float n,int k,float factor);
float gammp(float a, float x);
void gser(float *gamser, float a, float x, float *gln);
void gcf(float *gammcf, float a, float x, float *gln);
float gammln(float xx);




int main(void)
{


     srand48(time(NULL));
     srand(time(NULL));

clock_t Hinicio = clock(); 

//Indicamos el tipo de alfabeto=2 para el cual queremos hallar el índice

  int dimension;
  dimension=2;
  

//definimos dimension:numero alfabeto(2,4), valor_JS el valor JS max encontrado y size=tamaño_genoma
//devuelve significacion
  double significacion;
  double valor_JS,JSFF;
  int size;
  
 //  significacion=sl_m_d(dimension,valor_JS,size);



    //Abrimos el fichero y lo introducimos en un vector de memoria dinámica S[].
    //S[] es de tamaño size_genome y solo tiene entradas ACGT o N en caso de nucleótido desconocido 
    //Variables de entrada: fichero



 int size_genome;
 int frase1,frase_sig;
 int genes_introducidos=0;




 frase1=1;
 char *S;
 S=(char *)malloc(frase1*sizeof(char*));

 ifstream ficheroEntrada;
 string frase;

 ficheroEntrada.open ("Enterococcus_faecium.fa");
 frase_sig=0;
 bool region_genoma=false;
 int fN=0;
 bool empiezaN=true;
 int cuenta=0;
 int marcador=0;
 int N_final=0;
 int nuevo_c=0;
 int sumatorio;
 float rannd;
 char guardamos_x;
 int contenido=1;

 int size_c=1;
 int *cortes;
 cortes=(int *)malloc(size_c*sizeof(int*));
   int nuevo_corte=1;
   int corte_anterior=0;
   int num_cortes_cero=0;

 int aux_contador=0; //borrar


   int fA,fC,fG,fT,ffa,ffc,ffg,fft;
   fA=fG=fC=fT=0;
   ffa=ffc=fft=ffg=0;

 while(!ficheroEntrada.eof())
 {
    getline(ficheroEntrada, frase);
    frase_sig=frase_sig+frase.length();
   if(frase[0]=='>' || (frase[0]=='S' && frase[1]=='Q') || (frase[0]=='O' && frase[1]=='R'&& frase[2]=='I' && frase[3]=='G' && frase[4]=='I' && frase[5]=='N'))
    {  
      getline(ficheroEntrada, frase);
      region_genoma=true;
      frase_sig=frase.length();
    }
    if(region_genoma==true)
    {
         for (int i=0; i < frase.length(); i++){frase.at(i)=toupper(frase.at(i));}
         for(int i=0; i<frase.length();i++)
         {
        aux_contador++;
            if(frase.at(i)=='A' || frase.at(i)=='T')
            {

                  S=(char *)realloc(S,contenido * sizeof(char*));
                  contenido++;
                  S[genes_introducidos]='A';      
                  genes_introducidos++;  
                  fA++;
                  ffa++;
            }
            if(frase.at(i)=='C' || frase.at(i)=='G')
            {
                  S=(char *)realloc(S,contenido * sizeof(char*));
                  contenido++;
                  S[genes_introducidos]='C';      
                  genes_introducidos++; 
                  fC++; 
                  ffc++;
            }
            if(frase.at(i)=='N' || frase.at(i)=='R' || frase.at(i)=='Y' || frase.at(i)=='S' || frase.at(i)=='W' || frase.at(i)=='K' || frase.at(i)=='M'|| frase.at(i)=='V' || frase.at(i)=='H' || frase.at(i)=='D'|| frase.at(i)=='B' || frase.at(i)=='X')
            {

                if(empiezaN==true)
                {               
                  nuevo_c=0;
                  empiezaN=false;
                }
                fN++;              
 
            }  
                if(empiezaN==false)
                {
                  nuevo_c++;
                } 

                if(fN!=(nuevo_c))
                {
                      N_final=nuevo_c-1;
                      empiezaN=true;
                      fN=0;
                      nuevo_c=0;               
                   if(N_final<=10)       
                   {
                         guardamos_x=S[genes_introducidos-1];
                         genes_introducidos=genes_introducidos-1;
                         for(int kk=0;kk<N_final;kk++)
                         {                
                              rannd=rand()%4;

                              if(rannd==0)
                              {
                                  S=(char *)realloc(S,contenido * sizeof(char*));
                                  contenido++;
                                  S[genes_introducidos]='A';
                                  genes_introducidos++;
                                  fA++;
                                  ffa++;
                              }
                              if(rannd==1)
                              {
                                  S=(char *)realloc(S,contenido * sizeof(char*));
                                  contenido++;
                                  S[genes_introducidos]='C';
                                  genes_introducidos++;
                                 fC++;
                                 ffc++;
                              }
                              if(rannd==2)
                              {
                                  S=(char *)realloc(S,contenido * sizeof(char*));
                                  contenido++;
                                  S[genes_introducidos]='A';
                                  genes_introducidos++;
                                 fG++;
                                 ffg++;
                              }
                              if(rannd==3)
                              {
                                  S=(char *)realloc(S,contenido * sizeof(char*));
                                  contenido++;
                                  S[genes_introducidos]='C';
                                  genes_introducidos++;
                                 fT++;
                                 fft++;
                               }

                          }
                                  S[genes_introducidos]=guardamos_x;
                                  genes_introducidos++;
                   }
                   if(N_final>10)       
                   {

                     sumatorio=ffa+ffc+ffg+fft-1;

                     if(S[genes_introducidos-1]==0){ffa=ffa-1;}
                     if(S[genes_introducidos-1]==1){ffc=ffc-1;}
                     if(S[genes_introducidos-1]==2){ffg=ffg-1;}
                     if(S[genes_introducidos-1]==3){fft=fft-1;}

                       if(sumatorio>0)
                       {
                          cout<<"cortes de ACGT desde  "<<marcador<<"  hasta  "<<genes_introducidos-1<<endl;

                          nuevo_corte=corte_anterior+7;
                          cortes=(int *)realloc(cortes,nuevo_corte * sizeof(int));
                          cortes[corte_anterior]=marcador;
                          cortes[corte_anterior+1]=genes_introducidos-1;
                          cortes[corte_anterior+2]=0;                 
                          cortes[corte_anterior+3]=ffa;                
                          cortes[corte_anterior+4]=ffc;                  
                          cortes[corte_anterior+5]=ffg;                    
                          cortes[corte_anterior+6]=fft;                 
                          corte_anterior=nuevo_corte;
                          num_cortes_cero++;


                          marcador=genes_introducidos-1;
                          ffa=ffc=fft=ffg=0;

                       }

                     if(S[genes_introducidos-1]==0){ffa=ffa+1;}
                     if(S[genes_introducidos-1]==1){ffc=ffc+1;}
                     if(S[genes_introducidos-1]==2){ffg=ffg+1;}
                     if(S[genes_introducidos-1]==3){fft=fft+1;}


                   }
                } 
               
        

         }		

    }

 }
 ficheroEntrada.close();
 size_genome=genes_introducidos;

  sumatorio=ffa+ffc+ffg+fft;
  if(sumatorio>0)
  {
     cout<<"cortes de ACGT desde  "<<marcador<<"  hasta  "<<genes_introducidos<<endl;
                          nuevo_corte=corte_anterior+7;
                          cortes=(int *)realloc(cortes,nuevo_corte * sizeof(int));
                          cortes[corte_anterior]=marcador;
                          cortes[corte_anterior+1]=genes_introducidos;
                          cortes[corte_anterior+2]=0;                 
                          cortes[corte_anterior+3]=ffa;                
                          cortes[corte_anterior+4]=ffc;                  
                          cortes[corte_anterior+5]=ffg;                    
                          cortes[corte_anterior+6]=fft;                 
                          corte_anterior=nuevo_corte;
                          num_cortes_cero++;
  }



 cout<<"Tamaño genoma:  "<<size_genome<<"  "<<aux_contador<<endl;


  int lg;
  lg=size_genome;
  int N_matriz_cortes=num_cortes_cero;




  // Vamos a realizar los cortes

cout<<endl;
cout<<" Significación impuesta - SCC - SCC*sig - Cortes hallados  "<<endl;
cout<<endl;

 bool es;
 int dominio_i,dominio_f,can_i,can_f,can_corte,l,l1,l2;
 int contador_cortes,corte;
 double  JS,can_JS,H,H1,H2,HFF,JSF,HCC,HCX,JSaux;
 double  pA,pG,pC,pT,ppa,ppc,ppg,ppt;
 int A,C,G,T;
double d1,d2,dx,dxx;
int num_r;
int A1,A2,C1,C2,G1,G2,T1,T2;
double  sig,sig_r;
bool mas_cortes;
int sig_corte;
double HC;
int lc,lcc;
int can_A,can_C,can_G,can_T;
int frec_A,frec_C,frec_G,frec_T;
int ccA,ccG,ccC,ccT;
int frec_abs;
int l_sum;
int no_region_N;
double SCCmax=0;
double SIGmax;
double SCCmax_aux=0;
bool condicion=true;
int inicial;
int cortes_c;
bool cut;
int num_ex;


         pA=fA; 
         if(pA!=0)
         {
           pA=pA/lg;
           pA=pA*log2(pA);
         }
         pC=fC;
         if(pC!=0)
         {
           pC=pC/lg;
           pC=pC*log2(pC);
         }
         pG=fG;
         if(pG!=0)
         {
           pG=pG/lg;
           pG=pG*log2(pG);
         }
         pT=fT;
         if(pT!=0)
         {
           pT=pT/lg;
           pT=pT*log2(pT);
         } 


         HFF=-pA-pC-pG-pT;
   cout<<"Entropía máxima:  "<<HFF<<endl;

       JSaux=HFF;
       for(int ij=0; ij<N_matriz_cortes; ij++)
       {
         num_ex=ij*7;
         lcc=cortes[num_ex+1]-cortes[num_ex];
         pA=cortes[num_ex+3]; 
         if(pA!=0)
         {
           pA=pA/lcc;
           pA=pA*log2(pA);
         }
         pC=cortes[num_ex+4];
         if(pC!=0)
         {
           pC=pC/lcc;
           pC=pC*log2(pC);
         }
         pG=cortes[num_ex+5];
         if(pG!=0)
         {
           pG=pG/lcc;
           pG=pG*log2(pG);
         }
         pT=cortes[num_ex+6];
         if(pT!=0)
         {
           pT=pT/lcc;
           pT=pT*log2(pT);
         } 
        HCX=-pA-pC-pG-pT;
        dxx=lcc;
        dxx=dxx/lg;
        JSaux=JSaux-(dxx)*HCX;
                   
       }

cout<<"Cortes primeros:  "<<JSaux<<endl;




 sig=0.99;

 for(int ss=0; ss<10; ss++)
 {
       JSFF=HFF;
       num_cortes_cero=0;  
       for(int ij=0; ij<N_matriz_cortes; ij++)
       {
              num_ex=ij*7;
              frec_abs=cortes[num_ex+3]+cortes[num_ex+4]+cortes[num_ex+5]+cortes[num_ex+6];
              cortes_c=cortes[num_ex+1]-cortes[num_ex];
              cortes[num_ex+2]=1;
              if(frec_abs!=0 && cortes_c>=8)
              {
                  cortes[num_ex+2]=0;
                 num_cortes_cero++;
              }    
              if(frec_abs!=0 && cortes_c<8)
              {
         lcc=cortes[num_ex+1]-cortes[num_ex];
         pA=cortes[num_ex+3]; 
         if(pA!=0)
         {
           pA=pA/lcc;
           pA=pA*log2(pA);
         }
         pC=cortes[num_ex+4];
         if(pC!=0)
         {
           pC=pC/lcc;
           pC=pC*log2(pC);
         }
         pG=cortes[num_ex+5];
         if(pG!=0)
         {
           pG=pG/lcc;
           pG=pG*log2(pG);
         }
         pT=cortes[num_ex+6];
         if(pT!=0)
         {
           pT=pT/lcc;
           pT=pT*log2(pT);
         } 
        HCX=-pA-pC-pG-pT;
        dxx=lcc;
        dxx=dxx/lg;
        JSFF=JSFF-(dxx)*HCX;
              }      
       }

       cut=true;
       if(num_cortes_cero==0)
       {
          cut=false;
       } 

 mas_cortes=true;
 do{

      es=true;
      for(int kk=0; kk<N_matriz_cortes;kk++)
      {
          num_ex=kk*7;
          if(cortes[num_ex+2]==0 && es==true)
          {

             contador_cortes=num_ex;
             dominio_i=cortes[num_ex];
             dominio_f=cortes[num_ex+1];
             ccA=cortes[((num_ex)+3)];
             ccC=cortes[((num_ex)+4)];
             ccG=cortes[((num_ex)+5)];
             ccT=cortes[((num_ex)+6)];
 
             es=false;
             break;
          }
      }
   
  
   JS=0;
   corte=dominio_i+4;
   l=dominio_f-dominio_i;
  
         pA=ccA; 
         if(pA!=0)
         {
           pA=pA/l;
           pA=pA*log2(pA);
         }
         pC=ccC;
         if(pC!=0)
         {
           pC=pC/l;
           pC=pC*log2(pC);
         }
         pG=ccG;
         if(pG!=0)
         {
           pG=pG/l;
           pG=pG*log2(pG);
         }
         pT=ccT;
         if(pT!=0)
         {
           pT=pT/l;
           pT=pT*log2(pT);
         } 


         H=-pA-pC-pG-pT;
//

        A=C=G=T=0;

        for(int i=dominio_i;i<corte;i++)
        {
            if(S[i]=='A') {A++;}
            if(S[i]=='C'){C++;}
            if(S[i]=='G'){G++;}
            if(S[i]=='T'){T++;}
        }
 
    A=ccA-A;
    C=ccC-C;
    G=ccG-G;
    T=ccT-T;
 
    for(int ly=0; ly<l-7;ly++)
    {
       if(corte!=dominio_i+4)
       {
         if(S[corte-1]=='A')
         {
           A=A-1;
         }
         if(S[corte-1]=='C')
         {
           C=C-1;
         }
         if(S[corte-1]=='G')
         {
           G=G-1;
         }
         if(S[corte-1]=='T')
         {
           T=T-1;
         }
       }

        l2=dominio_f-corte;

       pA=A; 
       if(pA!=0)
       {
         pA=pA/l2;
         pA=pA*log2(pA);
       }
        pC=C;
       if(pC!=0)
       {
         pC=pC/l2;
         pC=pC*log2(pC);
       }
       pG=G;
       if(pG!=0)
       {
         pG=pG/l2;
         pG=pG*log2(pG);
       }
       pT=T;
       if(pT!=0)
       {
         pT=pT/l2;
         pT=pT*log2(pT);
       } 

       H2=-pA-pC-pG-pT;


              l1=corte-dominio_i;

         pA=(ccA-A); 
         if(pA!=0)
         {
           pA=pA/l1;
           pA=pA*log2(pA);
         }
         pC=ccC-C;
         if(pC!=0)
         {
           pC=pC/l1;
           pC=pC*log2(pC);
         }
         pG=ccG-G;
         if(pG!=0)
         {
           pG=pG/l1;
           pG=pG*log2(pG);
         }
         pT=ccT-T;
         if(pT!=0)
         {
           pT=pT/l1;
           pT=pT*log2(pT);
         } 

         H1=-pA-pC-pG-pT;



        d1=l1;
        d1=d1/l;
        d2=l2;
        d2=d2/l;

        can_JS=H-(d1)*H1-(d2)*H2;
        can_JS=can_JS;

        if(can_JS>=JS)
        {
            JS=can_JS;
            can_corte=corte;
            can_A=A;
            can_C=C;
            can_G=G;
            can_T=T;

        }


        
        corte++;
    }


//

// aquí vemos si hacemos o no el corte en can_corte dependiendo si cumple la condición de la significatividad
//sig_r por lo tanto corresponde a la Probabilidad(JS(random)>=JS(G))

     valor_JS=JS;
     size=l;
     significacion=sl_m_d(dimension,valor_JS,size);

      if(significacion>=sig)
      {

//HACEMOS CORTE y escribimos ambas subcadenas que salen del corte en la matriz cortes_dna, y borramos la cadena anterior
        
                  num_cortes_cero++;

                  cortes[contador_cortes]=dominio_i;
                  cortes[contador_cortes+1]=can_corte;
                  cortes[contador_cortes+2]=0;
                  frec_A=cortes[contador_cortes+3]-can_A;
                  cortes[contador_cortes+3]=frec_A;
                  frec_C=cortes[contador_cortes+4]-can_C;
                  cortes[contador_cortes+4]=frec_C;
                  frec_G=cortes[contador_cortes+5]-can_G;
                  cortes[contador_cortes+5]=frec_G;
                  frec_T=cortes[contador_cortes+6]-can_T;
                  cortes[contador_cortes+6]=frec_T;


                  if(can_corte-dominio_i<8)
                  {

                      cortes[contador_cortes+2]=1;
                      num_cortes_cero=num_cortes_cero-1;
         lcc=cortes[contador_cortes+1]-cortes[contador_cortes];
         pA=cortes[contador_cortes+3]; 
         if(pA!=0)
         {
           pA=pA/lcc;
           pA=pA*log2(pA);
         }
         pC=cortes[contador_cortes+4]; 
         if(pC!=0)
         {
           pC=pC/lcc;
           pC=pC*log2(pC);
         }
         pG=cortes[contador_cortes+5]; 
         if(pG!=0)
         {
           pG=pG/lcc;
           pG=pG*log2(pG);
         }
         pT=cortes[contador_cortes+6]; 
         if(pT!=0)
         {
           pT=pT/lcc;
           pT=pT*log2(pT);
         } 
        HCX=-pA-pC-pG-pT;
        dxx=lcc;
        dxx=dxx/lg;
        JSFF=JSFF-(dxx)*HCX;
                  }

               
                                 nuevo_corte=corte_anterior+7;

                       cortes=(int *)realloc(cortes,nuevo_corte * sizeof(int));
                
                     cortes[corte_anterior]=can_corte;
                     cortes[corte_anterior+1]=dominio_f;
                     cortes[corte_anterior+2]=0;
                     cortes[corte_anterior+3]=can_A;
                     cortes[corte_anterior+4]=can_C;
                     cortes[corte_anterior+5]=can_G;
                     cortes[corte_anterior+6]=can_T;
                  
               
           
 
                  

                  if(dominio_f-can_corte<8)
                  {
                      cortes[corte_anterior+2]=1;
                      num_cortes_cero=num_cortes_cero-1;
         lcc=cortes[corte_anterior+1]-cortes[corte_anterior];
         pA=cortes[corte_anterior+3]; 
         if(pA!=0)
         {
           pA=pA/lcc;
           pA=pA*log2(pA);
         }
         pC=cortes[corte_anterior+4];
         if(pC!=0)
         {
           pC=pC/lcc;
           pC=pC*log2(pC);
         }
         pG=cortes[corte_anterior+5];
         if(pG!=0)
         {
           pG=pG/lcc;
           pG=pG*log2(pG);
         }
         pT=cortes[corte_anterior+6];
         if(pT!=0)
         {
           pT=pT/lcc;
           pT=pT*log2(pT);
         } 
        HCX=-pA-pC-pG-pT;
        dxx=lcc;
        dxx=dxx/lg;
        JSFF=JSFF-(dxx)*HCX;

                  }
                  corte_anterior=nuevo_corte;
                  N_matriz_cortes=nuevo_corte/7;



                 
      
      }
      if(significacion<sig)
      {

// NO HACEMOS CORTE

         cortes[contador_cortes+2]=1;
         num_cortes_cero=num_cortes_cero-1;
         lcc=cortes[contador_cortes+1]-cortes[contador_cortes];
         pA=cortes[contador_cortes+3]; 
         if(pA!=0)
         {
           pA=pA/lcc;
           pA=pA*log2(pA);
         }
         pC=cortes[contador_cortes+4];
         if(pC!=0)
         {
           pC=pC/lcc;
           pC=pC*log2(pC);
         }
         pG=cortes[contador_cortes+5];
         if(pG!=0)
         {
           pG=pG/lcc;
           pG=pG*log2(pG);
         }
         pT=cortes[contador_cortes+6];
         if(pT!=0)
         {
           pT=pT/lcc;
           pT=pT*log2(pT);
         } 
        HCX=-pA-pC-pG-pT;
        dxx=lcc;
        dxx=dxx/lg;
        JSFF=JSFF-(dxx)*HCX;

                
      }

     cut=true;
     if(num_cortes_cero==0)
     {
         cut=false;
     }

 


 }while(cut==true);

//
//      cout<<endl;
//      cout<<"Número de cortes:  "<<N_matriz_cortes<<endl;
//      for(int kkk=0; kkk<N_matriz_cortes;kkk++)
//      {
//          cout<<cortes[kkk][1]<<"  ";
//      }
//      cout<<endl;   
    
 //
 


//Hallamos JS teniendo en cuenta ya todos los cortes que se han realizado.

int suma_cero;





                
          
      
      JSFF=JSFF-JSaux;
      SCCmax_aux=JSFF*sig;
      if(SCCmax<=SCCmax_aux)
      {
         SCCmax=SCCmax_aux;
         SIGmax=sig;
      }
      cout<<endl;
      cout<<sig<<"  "<<JSFF<<"  "<<JSFF*sig<<"  "<<N_matriz_cortes<<endl;

//Finalmente ordeno la matriz de cortes


//
//     inicial=0;
   
//     if(cortes[0][0]!=0)
//     {
//      for(int kk=0; kk<N_matriz_cortes;kk++)
//      {
//         if(cortes[kk][0]==0)
//         {
//            inicial=kk;
//            break;
//         }
  //    }       
//     }

//      for(int kk=0; kk<N_matriz_cortes;kk++)
//      {
//         for(int w=0; w<7; w++)
//         {
//           aux_cortes[w]=cortes[kk][w];
//         }
//         for(int w=0; w<7; w++)
//         {
//           cortes[kk][w]=cortes[inicial][w];
//         }
//         for(int w=0; w<7; w++)
//         {
//           cortes[inicial][w]=aux_cortes[w];
//         }

//            for(int kkk=0; kkk<N_matriz_cortes;kkk++)
//            {
//              if(cortes[kkk][0]==cortes[kk][1])
//              {
//                  inicial=kkk;
//                  break;
//              }
//            }
        
//      } 



//

sig=sig-0.01;
}

sig=sig+0.02;
free(S);
free(cortes);

cout<<endl;
SCCmax=SCCmax;
cout<<"Significación para la cual SCC es máxima  -   SCCmáxima  "<<endl;
cout<<SIGmax<<"   "<<SCCmax<<endl;
cout<<endl;

cout<<"Tiempo de ejecución:  "<<endl;
clock_t end = clock() - Hinicio; 
end=end/CLOCKS_PER_SEC;
cout<<end<<"  segundos"<<endl;
end=end/60;
cout<<end<<"  minutos"<<endl;






}





float sl_m_d(int dim,float d, int n)
{

   double sig;
   double ne;
   float f=0.8;

   if(dim==2)
   {
      ne=2.96*log(n)-7.88;
      sig=gammp_n_f(2*d*log(2)*n,ne,dim,f);
   }
   if(dim==4)
   {
      ne=2.44*log(n)-6.15;
      sig=gammp_n_f(2*d*log(2)*n,ne,dim,f);
   } 

   return sig;
}

float gammp_n_f(float x, float n,int k,float factor)
{
   float fsig,inter;
   float def1,def2;

   def1=(k-1);
   def1=def1/2;
   def2=factor*x;
   def2=def2/2;


   inter=gammp(def1,def2);
   if(inter<0.1)
   {
       fsig=0;
   }
   else
   {
       fsig=pow(inter,n);
   }
    
   return fsig;
}




float gammp(float a, float x)
// Returns  the  incomplete  gamma  function P(a,x).
{
   void gcf(float *gammcf, float a, float x, float *gln);
   void gser(float *gamser, float a, float x, float *gln);
   float gamser,gammcf,gln;
   if (x < 0.0 || a <= 0.0) 
   {
     cout<<"error in gammp"<<endl;
   }
   if (x < (a+1.0)) 
   {
     // Use  the  series  representation.
     gser(&gamser,a,x,&gln);
     return gamser;
   }
   else 
   {
      //Use the  continued  fraction  representation
     gcf(&gammcf,a,x,&gln);
     return 1.0-gammcf;
     //and  take  its  complement.
   }

}

void gser(float *gamser, float a, float x, float *gln)
   //Returns the incomplete gamma function P(a,x) evaluated by its series representation as gamser. Also  returns lnΓ(a) as gln.
{
   float gammln(float xx);
   int n;
   float sum,del,ap;

   *gln=gammln(a);

   if (x <= 0.0) 
   {
      if (x < 0.0)
      { 
         cout<<"error in gser"<<endl;
      }
         *gamser=0.0;
         return;
      
   }
   else 
   {
      ap=a;
      del=sum=1.0/a;
      for (n=1;n<=ITMAX;n++) 
      {
         ++ap;
         del *= x/ap;
         sum += del;
         if (fabs(del) < fabs(sum)*EPS) 
         {
              *gamser=sum*exp(-x+a*log(x)-(*gln));
              return;
         }
      }
      cout<<"error in gser else"<<endl;
      return;
   }
}












   //Returns  the  incomplete  gamma  function Q(a,x) evaluated  by  its  continued  fraction  representation  as gammcf.Also  returns 
   //lnΓ  (a)asgln.

void gcf(float *gammcf, float a, float x, float *gln)
{
   float gammln(float xx);
   int i;
   float an,b,c,d,del,h;

   *gln=gammln(a);
   b=x+1.0-a;

   //Set up for evaluating continued  fraction by  modified  Lentz’s  method with b=0.

   c=1.0/FPMIN;
   d=1.0/b;
   h=d;

   for (i=1;i<=ITMAX;i++) 
   {

     an = -i*(i-a);
     b += 2.0;
     d=an*d+b;
     if (fabs(d) < FPMIN) d=FPMIN;
     c=b+an/c;
     if (fabs(c) < FPMIN) c=FPMIN;
     d=1.0/d;
     del=d*c;
     h *= del;
     if (fabs(del-1.0) < EPS) break;
   }
   if (i > ITMAX) 
   {
      cout<<"error a too large, ITMAX too small in gcf"<<endl;
   }
   *gammcf=exp(-x+a*log(x)-(*gln))*h;

}


   //Returns the value ln[Γ(xx)] for xx>0.

float gammln(float xx)
{

   double x,y,tmp,ser;
   static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
   int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}








