/**************************************************************************************/  
//  Computes the Sequence Compositional Complexity (SCC_SW) of a genome 
//  (https://link.aps.org/doi/10.1103/PhysRevLett.80.1344)
//  Related works:  -> https://doi.org/10.3390/biology12020322 
//                  -> https://doi.org/10.1038/s41598-020-76014-4
/**************************************************************************************/  
//  Code written by Rebeca de la Fuente :  2022                                 
/**************************************************************************************/  
//** Output Notes:
//** Este programa saca por pantalla:
//   Tamaño genoma
//Cortes correspondientes a los nucleótidos desconocidos,'N'
// Nivel de significación impuesta - SCC - SCC*sig - número de cortes          **\\
// Nivel de sig para la cual la SCC se maximiza  -  SCCmax


#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<cstdlib>
#include<string.h>
#include<fstream>
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

//Indicamos el tipo de alfabeto para el cual queremos hallar el índice

  int dimension;
  dimension=2;
  

//definimos dimension:numero alfabeto(2,4), valor_JS el valor JS max encontrado y size=tamaño_genoma
//devuelve significacion
  double significacion;
  double valor_JS;
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

 ficheroEntrada.open ("prueba4.adn");
 frase_sig=0;
 bool region_genoma=false;



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
         S=(char *)realloc(S,frase_sig * 2 * sizeof(char*));
         for (int i=0; i < frase.length(); i++){frase.at(i)=toupper(frase.at(i));}
         for(int i=0; i<frase.length();i++)
         {
            if(frase.at(i)=='C'  || frase.at(i)=='G')
            {
                  S[genes_introducidos]='C';      
                  genes_introducidos++;  
            }
            if(frase.at(i)=='A' || frase.at(i)=='T')
            {
                  S[genes_introducidos]='A';      
                  genes_introducidos++;  
            }
            if(frase.at(i)=='N' || frase.at(i)=='R' || frase.at(i)=='Y' || frase.at(i)=='S' || frase.at(i)=='W' || frase.at(i)=='K' || frase.at(i)=='M'|| frase.at(i)=='V' || frase.at(i)=='H' || frase.at(i)=='D'|| frase.at(i)=='B' || frase.at(i)=='X')
            {
                  S[genes_introducidos]='N';      
                  genes_introducidos++;  
            }         

         }		

    }

 }
 ficheroEntrada.close();
 S=(char *)realloc(S,genes_introducidos * 2 * sizeof(char*));
 size_genome=genes_introducidos;

 cout<<"Tamaño genoma:  "<<size_genome<<endl;

    //Definimos la matriz de cortes[][] a través de memoria dinámica. Cada fila corresponde a una subcadena 
    //cortes[][0] define posición inicial al corte
    //cortes[][1] define posición final+1 al corte
    //cortes[][2] define '0' si falta por visitar y '1' si ya no hay más cortes
    //cortes[][3-6] las frecuencias absolutas de A,C,G,T

   int **cortes;
   int size_c=7;
   int row1=1;

   cortes=(int **)malloc(row1*sizeof(int*));
   for(int i=0; i<row1; i++)
   {
         cortes[i]=(int*)malloc((size_c)*sizeof(int));
   }

   //Variables de frecuencias absolutas
  
   int fA,fC,fG,fT;

   //Arreglo la cadena con respecto a los carácteres N y preparo cortes de entrada



   int fN;
   int N_inicio,N_final;
   int nuevo_c;
   fA=fG=fC=fT=0;
   int recorremos_S=0;
   int r_anterior=0;
   int nuevo_corte=1;
   int corte_anterior=0;
   float rannd;
   do{
        
        if(S[recorremos_S]=='N')
        {
           fN=0;
           N_inicio=recorremos_S;
           nuevo_c=recorremos_S;
           do{
              if(S[nuevo_c]=='N')
              {
                 fN++;
              }
              nuevo_c++;
             }while(nuevo_c<=size_genome && fN==(nuevo_c-N_inicio));
           N_final=nuevo_c-1;
           recorremos_S=N_final;  
      

          if(N_final-N_inicio<=10)       
          {
             for(int kk=N_inicio;kk<N_final;kk++)
             {
               rannd=rand()%4;

               if(rannd==0)
               {
                 S[kk]='C';
                 fA++;
               }
               if(rannd==1)
               {
                 S[kk]='A';
                 fC++;
               }
               if(rannd==2)
               {
                 S[kk]='C';
                 fG++;
               }
               if(rannd==3)
               {
                S[kk]='A';
                fT++;
               }
             }
          }
          if(N_final-N_inicio>10)       
          {
             r_anterior=(fA+fC+fG+fT);
             if((r_anterior)==0)
             {
                cout<<"cortes de N desde: "<<N_inicio<<"  hasta  "<<N_final<<endl;
                     
                      nuevo_corte=corte_anterior+1;
                      cortes = (int**) realloc(cortes, nuevo_corte * 2 * sizeof (int));
                      for(int iii=corte_anterior; iii<nuevo_corte; iii++)
                      {
                          cortes[iii]=(int*)malloc((size_c)*sizeof(int));
                      }
                
                     cortes[corte_anterior][0]=N_inicio;
                     cortes[corte_anterior][1]=N_final;
                     cortes[corte_anterior][2]=1;
                     cortes[corte_anterior][3]=0;
                     cortes[corte_anterior][4]=0;
                     cortes[corte_anterior][5]=0;
                     cortes[corte_anterior][6]=0;
                     corte_anterior=nuevo_corte;
             }

             if((r_anterior)>0)
             {
                  cout<<"cortes de ACGT desde  "<<(N_inicio-r_anterior)<<"  hasta  "<<N_inicio<<endl;
                      nuevo_corte=corte_anterior+1;
                      cortes = (int**) realloc(cortes, nuevo_corte * 2 * sizeof (int));
                      for(int iii=corte_anterior; iii<nuevo_corte; iii++)
                      {
                          cortes[iii]=(int*)malloc((size_c)*sizeof(int));
                      }
                
                     cortes[corte_anterior][0]=(N_inicio-r_anterior);
                     cortes[corte_anterior][1]=N_inicio;
                     cortes[corte_anterior][2]=0;
                     cortes[corte_anterior][3]=fA;
                     cortes[corte_anterior][4]=fC;
                     cortes[corte_anterior][5]=fG;
                     cortes[corte_anterior][6]=fT;
                     corte_anterior=nuevo_corte;
                     fA=fC=fG=fT=0;
                  cout<<"cortes de N desde  "<<(N_inicio)<<"  hasta  "<<N_final<<endl;
                      nuevo_corte=corte_anterior+1;
                      cortes = (int**) realloc(cortes, nuevo_corte * 2 * sizeof (int));
                      for(int iii=corte_anterior; iii<nuevo_corte; iii++)
                      {
                          cortes[iii]=(int*)malloc((size_c)*sizeof(int));
                      }
                
                     cortes[corte_anterior][0]=N_inicio;
                     cortes[corte_anterior][1]=N_final;
                     cortes[corte_anterior][2]=1;
                     cortes[corte_anterior][3]=0;
                     cortes[corte_anterior][4]=0;
                     cortes[corte_anterior][5]=0;
                     cortes[corte_anterior][6]=0;
                     corte_anterior=nuevo_corte;
             }
          }

                      
        }

 
       
         if(S[recorremos_S]=='A'){fA++;}
         if(S[recorremos_S]=='C'){fC++;}
         if(S[recorremos_S]=='G'){fG++;}
         if(S[recorremos_S]=='T'){fT++;}  

        recorremos_S++;
     }while(recorremos_S<size_genome);

      
             r_anterior=(fA+fC+fG+fT);
             if((r_anterior)>0)
             {
                  cout<<"cortes de ACGT desde  "<<(size_genome-r_anterior)<<"  hasta  "<<size_genome<<endl;
                      nuevo_corte=corte_anterior+1;
                      cortes = (int**) realloc(cortes, nuevo_corte * 2 * sizeof (int));
                      for(int iii=corte_anterior; iii<nuevo_corte; iii++)
                      {
                          cortes[iii]=(int*)malloc((size_c)*sizeof(int));
                      }
                
                     cortes[corte_anterior][0]=(size_genome-r_anterior);
                     cortes[corte_anterior][1]=size_genome;
                     cortes[corte_anterior][2]=0;
                     cortes[corte_anterior][3]=fA;
                     cortes[corte_anterior][4]=fC;
                     cortes[corte_anterior][5]=fG;
                     cortes[corte_anterior][6]=fT;
                     corte_anterior=nuevo_corte;
                  fA=fC=fG=fT=0;
             }
         


  fA=fC=fG=fT=0;


  int N_matriz_cortes=nuevo_corte;





  // Vamos a realizar los cortes

cout<<endl;
cout<<" Significación impuesta - SCC - SCC*sig - Cortes hallados  "<<endl;
cout<<endl;


 bool es;
 int dominio_i,dominio_f,can_i,can_f,can_corte,l,l1,l2;
 int contador_cortes,corte;
 double JS,can_JS,H,H1,H2;
 double pA,pG,pC,pT;
 int A,C,G,T;
double d1,d2;
int num_r;
int A1,A2,C1,C2,G1,G2,T1,T2;
float sig,sig_r;
bool mas_cortes;
int sig_corte;
double HC;
int lc;
int can_A,can_C,can_G,can_T;
int frec_A,frec_C,frec_G,frec_T;
int ccA,ccG,ccC,ccT;
int frec_abs;
int l_sum;
int no_region_N;
float SCCmax=0;
float SIGmax;
float SCCmax_aux=0;
bool condicion=true;
int inicial;
int aux_cortes[7];


 sig=0.99;

 for(int ss=0; ss<49; ss++)
 {

       for(int ij=0; ij<N_matriz_cortes; ij++)
       {
              frec_abs=cortes[ij][3]+cortes[ij][4]+cortes[ij][5]+cortes[ij][6];
              if(frec_abs!=0)
              {
                 cortes[ij][2]=0;
              }          
       } 

 mas_cortes=true;
 do{

      es=true;
      for(int kk=0; kk<N_matriz_cortes;kk++)
      {
          if(cortes[kk][2]==0 && es==true)
          {
             contador_cortes=kk;
             dominio_i=cortes[kk][0];
             dominio_f=cortes[kk][1];
             ccA=cortes[kk][3];
             ccC=cortes[kk][4];
             ccG=cortes[kk][5];
             ccT=cortes[kk][6];
 
             es=false;
          }
      }

      condicion=false;
      if((dominio_f-dominio_i)>=4)
      {
          condicion=true;
      }

  if(condicion==false)
  {
      cortes[contador_cortes][2]=1;
  }

   if(condicion==true)
   {
   
   

   JS=0;
   corte=dominio_i+1;
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
 
    for(int ly=0; ly<l-1;ly++)
    {

        l2=dominio_f-corte;
        A=C=G=T=0;
        for(int i=corte;i<dominio_f;i++)
        {
            if(S[i]=='A') {A++;}
            if(S[i]=='C'){C++;}
            if(S[i]=='G'){G++;}
            if(S[i]=='T'){T++;}
        }
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

        if(can_JS>=JS && (corte-4)>=dominio_i && (corte+4)<=dominio_f)
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
        
         sig_corte=(can_corte-4);
         if(sig_corte!=dominio_i)
         {           
                  cortes[contador_cortes][0]=dominio_i;
                  cortes[contador_cortes][1]=can_corte;
                  cortes[contador_cortes][2]=0;
                  frec_A=cortes[contador_cortes][3]-can_A;
                  cortes[contador_cortes][3]=frec_A;
                  frec_C=cortes[contador_cortes][4]-can_C;
                  cortes[contador_cortes][4]=frec_C;
                  frec_G=cortes[contador_cortes][5]-can_G;
                  cortes[contador_cortes][5]=frec_G;
                  frec_T=cortes[contador_cortes][6]-can_T;
                  cortes[contador_cortes][6]=frec_T;

                  
               
           
         }
         if(sig_corte==dominio_i)
         {
                  cortes[contador_cortes][0]=dominio_i;
                  cortes[contador_cortes][1]=can_corte;
                  cortes[contador_cortes][2]=1;
                  frec_A=cortes[contador_cortes][3]-can_A;
                  cortes[contador_cortes][3]=frec_A;
                  frec_C=cortes[contador_cortes][4]-can_C;
                  cortes[contador_cortes][4]=frec_C;
                  frec_G=cortes[contador_cortes][5]-can_G;
                  cortes[contador_cortes][5]=frec_G;
                  frec_T=cortes[contador_cortes][6]-can_T;
                  cortes[contador_cortes][6]=frec_T;


                  
         }


         sig_corte=(can_corte+4);
         if(sig_corte!=dominio_f)
         {
                      nuevo_corte=corte_anterior+1;
                      cortes = (int**) realloc(cortes, nuevo_corte * 2 * sizeof (int));
                      for(int iii=corte_anterior; iii<nuevo_corte; iii++)
                      {
                          cortes[iii]=(int*)malloc((size_c)*sizeof(int));
                      }               
                  cortes[corte_anterior][0]=can_corte;
                  cortes[corte_anterior][1]=dominio_f;
                  cortes[corte_anterior][2]=0;
                  cortes[corte_anterior][3]=can_A;
                  cortes[corte_anterior][4]=can_C;
                  cortes[corte_anterior][5]=can_G;
                  cortes[corte_anterior][6]=can_T;
                  corte_anterior=nuevo_corte;
                  N_matriz_cortes=nuevo_corte;

                 
         }
         if(sig_corte==dominio_f)
         {
                      nuevo_corte=corte_anterior+1;
                      cortes = (int**) realloc(cortes, nuevo_corte * 2 * sizeof (int));
                      for(int iii=corte_anterior; iii<nuevo_corte; iii++)
                      {
                          cortes[iii]=(int*)malloc((size_c)*sizeof(int));
                      }               
                  cortes[corte_anterior][0]=can_corte;
                  cortes[corte_anterior][1]=dominio_f;
                  cortes[corte_anterior][2]=1;
                  cortes[corte_anterior][3]=can_A;
                  cortes[corte_anterior][4]=can_C;
                  cortes[corte_anterior][5]=can_G;
                  cortes[corte_anterior][6]=can_T;
                  corte_anterior=nuevo_corte;
                  N_matriz_cortes=nuevo_corte;

                  
         }
         
      }
      if(significacion<sig)
      {

// NO HACEMOS CORTE

         cortes[contador_cortes][2]=1;

                
      }



      mas_cortes=false;
      for(int kk=0; kk<N_matriz_cortes;kk++)
      {
          if(cortes[kk][2]==0)
          {
             mas_cortes=true;
          }
      }

  }


 }while(mas_cortes==true);

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

       l_sum=0;
       A=C=G=T=0;
       for(int i=0; i<N_matriz_cortes; i++)
       {
          A=A+cortes[i][3];
          C=C+cortes[i][4];
          G=G+cortes[i][5];
          T=T+cortes[i][6];
          if(cortes[i][3]==0 && cortes[i][4]==0 && cortes[i][5]==0 && cortes[i][6]==0)
          {
             l_sum=l_sum+(cortes[i][1]-cortes[i][0]);
          }          
       } 
       l=size_genome-l_sum;

         pA=A; 
         if(pA!=0)
         {
           pA=pA/l;
           pA=pA*log2(pA);
         }
         pC=C;
         if(pC!=0)
         {
           pC=pC/l;
           pC=pC*log2(pC);
         }
         pG=G;
         if(pG!=0)
         {
           pG=pG/l;
           pG=pG*log2(pG);
         }
         pT=T;
         if(pT!=0)
         {
           pT=pT/l;
           pT=pT*log2(pT);
         } 


         H=-pA-pC-pG-pT;
         JS=H;


      for(int kk=0; kk<N_matriz_cortes;kk++)
      {
       no_region_N=cortes[kk][3]+cortes[kk][4]+cortes[kk][5]+cortes[kk][6];
       if(no_region_N!=0)
       {
               dominio_i=cortes[kk][0];
               dominio_f=cortes[kk][1];
               lc=dominio_f-dominio_i;
               A=cortes[kk][3];
               C=cortes[kk][4];
               G=cortes[kk][5];
               T=cortes[kk][6];
        
         pA=A; 
         if(pA!=0)
         {
           pA=pA/lc;
           pA=pA*log2(pA);
         }
         pC=C;
         if(pC!=0)
         {
           pC=pC/lc;
           pC=pC*log2(pC);
         }
         pG=G;
         if(pG!=0)
         {
           pG=pG/lc;
           pG=pG*log2(pG);
         }
         pT=T;
         if(pT!=0)
         {
           pT=pT/lc;
           pT=pT*log2(pT);
         } 


         HC=-pA-pC-pG-pT;

        d1=lc;
        d1=d1/l;
    

        JS=JS-(d1)*HC;
      }          
          
      }

      SCCmax_aux=JS*sig;
      if(SCCmax<=SCCmax_aux)
      {
         SCCmax=SCCmax_aux;
         SIGmax=sig;
      }
      cout<<endl;
      cout<<sig<<"  "<<JS<<"  "<<JS*sig<<"  "<<N_matriz_cortes<<endl;

//Finalmente ordeno la matriz de cortes



     inicial=0;
   
     if(cortes[0][0]!=0)
     {
      for(int kk=0; kk<N_matriz_cortes;kk++)
      {
         if(cortes[kk][0]==0)
         {
            inicial=kk;
            break;
         }
      }       
     }

      for(int kk=0; kk<N_matriz_cortes;kk++)
      {
         for(int w=0; w<7; w++)
         {
           aux_cortes[w]=cortes[kk][w];
         }
         for(int w=0; w<7; w++)
         {
           cortes[kk][w]=cortes[inicial][w];
         }
         for(int w=0; w<7; w++)
         {
           cortes[inicial][w]=aux_cortes[w];
         }

            for(int kkk=0; kkk<N_matriz_cortes;kkk++)
            {
              if(cortes[kkk][0]==cortes[kk][1])
              {
                  inicial=kkk;
                  break;
              }
            }
        
      } 




//


sig=sig-0.02;
}

sig=sig+0.02;

cout<<endl;
cout<<"Significación para la cual SCC es máxima  -   SCCmáxima  "<<endl;
cout<<SIGmax<<"   "<<SCCmax<<endl;
cout<<endl;

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










