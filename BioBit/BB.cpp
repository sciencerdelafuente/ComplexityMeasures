/**************************************************************************************/  
//  Computes the BioBit (/doi.org/10.1038/srep28840) from a genome S[L] of size L.
//  Related works:  -> https://doi.org/10.3390/biology12020322 
//                 -> https://doi.org/10.1038/s41598-020-76014-4
/**************************************************************************************/  
//  Code written by Rebeca de la Fuente :  2020                                 
/**************************************************************************************/ 
#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<cstdlib>

using namespace std;

int main(void)
{


  srand48(time(NULL));
  srand(time(NULL));

// L es el tamaño 'n' del genoma y S[L] el genoma modelo del programa.
// Defino LG2 como log2(L) y hallo los enteros k1 y k2 entre los q se encuentra.

int L=13;
char S[L]={'A','C','T','A','C','T','A','G','C','C','A','C','T'};

float LG2;
int k1,k2;
int dic1,dic2;
int contador,contador_k;
bool esta;


LG2=L;
LG2=log2(LG2);
k1=(int)LG2;

LG2=L;
LG2=log2(LG2)+1;
k2=(int)LG2;

LG2=L;
LG2=log2(LG2);

dic1=L-k1+1;
dic2=L-k2+1;

// Los vectores num, cada uno para su k correspondiente, darán lugar al final del programa al número de
// ocurrencias de cada una de las palabras existentes en el genoma, de forma que hay un total de dic=n-k+1
// número máximo de palabras que existen en el genoma.

int num1[dic1];
int num2[dic2];

 for(int i=0;i<dic1;i++)
 {
    num1[i]=1;
 }
 for(int i=0;i<dic2;i++)
 {
    num2[i]=1;
 }

char v1[dic1][k1];
char v2[dic2][k2];


//Almaceno en v[][] todas las palabras del genoma, de forma que cada fila
// está compuesta por una palabra del genoma. Me dan igual las repeticiones, poque más tarde lo único que haré será 
// borrar las repetidas y apuntar las ocurrencias de cada una de ellas en num[]. 
// la posición x del vector num[x], corresponderá a la palabra que se encuentra en la fila x de la matriz v[][]

 contador=0;

 for(int i=0;i<dic1;i++)
 {

     for(int j=0; j<k1;j++)
     {
         v1[i][j]=S[contador];       
         contador++;

     }
     contador=contador-k1+1;
 }

 contador=0;

 for(int i=0;i<dic2;i++)
 {

     for(int j=0; j<k2;j++)
     {
         v2[i][j]=S[contador];       
         contador++;

     }
     contador=contador-k2+1;
 }

//OBTENEMOS v[][] como la matriz formada por todas las palabras sin repeticiones, y num[] el vector que contendrá en num. de ocurrencias
//de cada palabra.
// Hay un total de dic1 palabras, que van de la fila 0 a la fila dic1-1
// Voy a ir escogiendo desde la palabra correspondiente a la fila 0 a la dic1-2 (primer bucle)
// , y cada una de ellas la compararé con las sucesivas
// hasta llegar a dic1-1 (segundo bucle). Esto lo hago para ver si hay palabras repetidas. SI hay alguna repetida, la borro de la matriz,
// reemplazando las letras que componen la palabra por 'N', y añadiré uno al número de ocurrencias de dicha palabra. 
// Por lo tanto 'contador' hará referencia a la palabra que compararé con todas las demás (las sucesivas hasta el final)
// y 'contador_k' serán todas las palabras con las que comparo la palabra correspondiente a 'contador'.

 for(int i=0;i<dic1-1;i++)
 {
  contador=i;
  contador_k=contador+1;
  for(int j=contador_k;j<dic1;j++)
  {

   if(v1[contador][0]!='N' || v1[j][0]!='N')
   {

      esta=true;
      for(int kk=0;kk<k1;kk++)
      {
         if(v1[contador][kk]!=v1[j][kk])
         {
            esta=false;
         }
      }

      if(esta==true)
      {
        for(int kk=0;kk<k1;kk++)
        {
          v1[j][kk]='N';
        } 
          num1[contador]++;
          num1[j]=0;
      }

    }
          
  } 
   
 }







 for(int i=0;i<dic2-1;i++)
 {
  contador=i;
  contador_k=contador+1;
  for(int j=contador_k;j<dic2;j++)
  {

   if(v2[contador][0]!='N' || v2[j][0]!='N')
   {

      esta=true;
      for(int kk=0;kk<k2;kk++)
      {
         if(v2[contador][kk]!=v2[j][kk])
         {
            esta=false;
         }
      }

      if(esta==true)
      {
        for(int kk=0;kk<k2;kk++)
        {
          v2[j][kk]='N';
        } 
          num2[contador]++;
          num2[j]=0;
      }

    }
          
  } 
   
 }

//Ahora que ya tenemos el número de ocurrencias, hallamos E1 y 	E2

float E1,E2,E;
float p;



  E1=0;
  for(int i=0;i<dic1;i++)
  {
    if(num1[i]!=0)
    {
       p=num1[i];
       p=p/dic1;
       p=p*log2(p);
       E1=E1-p;
    }
  }

  E2=0;
  for(int i=0;i<dic2;i++)
  {
    if(num2[i]!=0)
    {
       p=num2[i];
       p=p/dic2;
       p=p*log2(p);
       E2=E2-p;
    }
  }

//Interpolamos los valores de (k1,E1) y (k2,E2) para hallar E sobre k=log2(n)

float LG;
LG=L;
LG=log(LG)/log(4);


E=(E2-E1);
E=E/(k2-k1);
E=E*(LG2-k1);
E=E+E1;

float AF,AC,BB;


//Hallamos BB(G)

AC=LG2-E;
AF=AC/LG;
BB=sqrt(LG)*sqrt(AF)*pow((1-2*AF),3);



cout<<BB<<endl;















}
