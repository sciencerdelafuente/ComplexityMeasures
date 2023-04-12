/**************************************************************************************/  
//  Computes the Genome Signature from a genome in FASTA format.
//  -> k-mer library 
//  -> Coarse-graining and k-mer frequencies into heat map
//  -> Genome Signature
//  Related works:  -> https://doi.org/10.3390/biology12020322 
//                  -> https://doi.org/10.1038/s41598-020-76014-4
/**************************************************************************************/  
//  Code written by Rebeca de la Fuente :  2022                                 
/**************************************************************************************/  
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iterator>
#include <map>
#include <sys/stat.h>
  
// g++ -std=c++11 genomeSig.cpp -o o
using namespace std;

  struct  position{
  double x;
  double y;
  };


void lectura(vector<char> &,string);
map<string,int> frecuencias_kmeros(vector<char> &,int);
void GenomeSignature(vector<char> &,int, int,string &);
void Mapping(vector<char> &,int &,string &);
void VertexConversion(char &,position &);


int main()
{

    string name="homosapiens"; 
    string arch=name+".fna";

    vector<char> genome;
    lectura(genome,arch);

    string stringpath = "Output_"+name; 
    mkdir(stringpath.c_str(),0777);


    int k=9;
    int kmin,kmax;
    kmin=1;
    kmax=8;
    //GenomeSignature(genome,kmin,kmax,stringpath);


    Mapping(genome,k,stringpath);


}




void Mapping(vector <char> &genome,int &k,string &stringpath)
{

    //string path="/home/rebeca/Escritorio/GS/"+stringpath+"/PointDensity.txt";
    //ofstream Outfile(path);
    string path="/home/rebeca/Escritorio/GS/"+stringpath+"/kDensity.txt";
    ofstream Outfile2(path);

    map<int,position> Coordinates;
    position ant, sig, vertex;
    ant.x=0.5;
    ant.y=0.5;




    int boxes=pow(2,k);
    double boxd=boxes;
    boxes=boxes+1;
    int L[boxes][boxes];
    double h;
    h=1/boxd;
    int boxx,boxy;

    for(int i=0;i<boxes;i++)
    {
    for(int j=0;j<boxes;j++)
    {
        L[i][j]=0;
    }
    }

    double dx,dy;
    

    for(int i=0;i<genome.size();i++)
    {
        VertexConversion(genome[i],vertex);
        sig.x=(ant.x+vertex.x)*0.5;
        sig.y=(ant.y+vertex.y)*0.5;
        //Coordinates.insert({i,sig});


        dx=sig.x;
        dy=sig.y;

        dx=dx*(boxd);
        boxx=dx;
        dy=dy*(boxd);
        boxy=dy;

        L[boxx][boxy]++;

        ant.x=sig.x;
        ant.y=sig.y;
        
    }

    // k-AVERAGE
    for(int i=0;i<boxes;i++)
    {
    for(int j=0;j<boxes;j++)
    {
        Outfile2<<(i+0.5)*h<<" "<<(j+0.5)*h<<" "<<L[i][j]<<endl;
    }
    }

    /*
    for(auto it = Coordinates.cbegin(); it != Coordinates.cend(); ++it)
    {
     //   Outfile<<it->first<<" "<< it->second.x <<" "<<it->second.y<<endl;
    }
    */



    // POINT DENSITY
    /*
    for(auto it = Coordinates.cbegin(); it != Coordinates.cend(); ++it)
    {
        dx=it->second.x;
        dy=it->second.y;

        dx=dx*(boxd);
        boxx=dx;
        dy=dy*(boxd);
        boxy=dy;

        Outfile<<it->second.x<<" "<<it->second.y<<"   "<<L[boxx][boxy]<<endl;
    }
    */
    



}





























void VertexConversion(char &ch,position &pos)
{
   
     if(ch=='A')
     {
        pos.x=0;
        pos.y=0;
     }

     if(ch=='C')
     {
        pos.x=0;
        pos.y=1;
     }

     if(ch=='G')
     {
        pos.x=1;
        pos.y=1;
     }

     if(ch=='T')
     {
        pos.x=1;
        pos.y=0;
     }

}







void GenomeSignature(vector<char> &genome,int kmin, int kmax, string &stringpath)
{


    map<int,map<string,int>> library;
    int G=genome.size();
    int N,n,number_kmers;
    int f;
    double expected_f;
    double GenomicSig;
    double GenomeSigs[kmax-kmin];


    string path="/home/rebeca/Escritorio/GS/"+stringpath+"/Genomic_Signature.txt";
    ofstream Outfile(path);

   

    Outfile<<"********* Genomic Signature **********"<<endl;
    Outfile<<endl;
    Outfile<<"Genome Size (G): "<<G<<endl;
    Outfile<<endl;
    int i=0;
    for(int k=kmin;k<=kmax;k++)
    {

        library[k]=frecuencias_kmeros(genome,k); 
        N=pow(4,k);   
        n=library[k].size();
        number_kmers=G-k+1;
        expected_f=number_kmers/N;
        //Genomic Signature ok k-mers
        GenomicSig=0;
        for(auto it = library[k].cbegin(); it != library[k].cend(); ++it)
        {
              f=it->second;
              GenomicSig=GenomicSig+fabs(1-(f/expected_f));
        }
        GenomicSig=GenomicSig+(N-n);
        GenomicSig=GenomicSig/N;
        GenomeSigs[i]=GenomicSig;
        i++;
       

        Outfile <<"k:"<<k<<endl;
        Outfile<<"   Total possible k-mers (N=4^k):"<<N<<endl;
        Outfile<<"   Num. of k-mers in Genome (n=G-k+1): "<<number_kmers<<endl;
        Outfile<<"   k-mers in Genome:"<<n<<endl;
        Outfile<<endl;
        Outfile<<"Genomic Signature: "<<GenomicSig<<endl;
        Outfile<<endl;
        Outfile<<endl;
        Outfile<<"Frequencies of k-mer Library"<<endl;
        Outfile<<endl;
        Outfile<<endl;
        for(auto it = library[k].cbegin(); it != library[k].cend(); ++it)
        {
              Outfile <<"    "<< it->first << " " << it->second <<endl;
        }
        Outfile<<endl;
        Outfile<<endl;
    }

}


map<string,int> frecuencias_kmeros(vector<char> &genome,int k)
{

    map<string,int> library;

    string s;  //Initial K-mer
    string key_to_find;

    int number_kmers=genome.size()-k+1;


    for(int j=0;j<k;j++)
    {
        s=s+genome[j];
    }

    for(int i=0;i<number_kmers;i++)
    {   
        //cout<<s<<" "; 
        key_to_find=s;
        if(library.find(key_to_find) != library.end())
        {
           //cout<<"esta"<<endl;
           library.find(key_to_find)->second++;
        }
        else
        {  
            //cout<<"no esta"<<endl;
            library.insert(pair<string, int>(key_to_find, 1));
        }
        s.erase(0,1);
        s=s+genome[i+k];
    }

    //for(auto it = library.cbegin(); it != library.cend(); ++it)
    //{
      // std::cout << it->first << " " << it->second <<endl;
    //}

    return library;


}





void lectura(vector<char> &lec,string name_arch)
{

    char w;
    string Initial;
    ifstream file(name_arch.c_str());  

    file >> w;
    if(w=='>')
    {
       getline(file,Initial);
       file >> w;
    }    

    do 
    {
       if(w=='A' || w=='C' || w=='G' || w=='T')
       {
          lec.push_back(w);
       }
       file >> w;

    }while(!file.eof());

    //for(int i=0;i<lec.size();i++)
    //{
    //    cout<<lec[i]<<" ";
    //}

    file.close();


}


