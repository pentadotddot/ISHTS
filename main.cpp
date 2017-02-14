#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <time.h>  /* time */
#include <math.h> //log iylesmi
#include <random>
#include <sys/types.h>
#include <unistd.h>
#include <sstream>
#include <chrono>
//#define _BSD_SOURCE
#include <sys/time.h>
//SZÜKSÉGES FÁJLOK:adatok.txt
//idő
//szintek száma
//mutrata (jelen verzióban mut várható száma)
//delta Nk pk

using namespace std;

struct history{
    int szint;
    string type;
    double idopont;
};

struct sejt
{
    int oszt_szam,mut_szam,n_acd,n_scd,n_scdplusd,id;
    double rscdplusd,rscd,racd;
    vector<history> historyVector;
};


struct szint //Ebben tárolom az egy szinthez tartozó rátákat és a szint várt méretét
{
    double rscdplusd,rscd,racd,delta,pk,qk;
    int Nk;
};

struct lista
{
    double ertek;
    int szintszam, sejtszam, tipus;
};



void txtsorcsere(string ezt, string erre, ifstream& in, const char *inStreamFileName);
void fancyRead(ifstream & ifs, string & tmp);

map<int,vector<sejt> > x;//ebbe tároljuk a sejteket
//map<string,vector<sejt> > x;//ebbe tároljuk a sejteket
//vector<double> idovektord;

vector<int> Histog_Counter_Varied;

//vector<double> idovektorD;121
//vector<int> Counter;


double idoskala;//egy idõlépés hossza
double ido; //amennyi ideig fut a szim.
double legyartando;
vector<szint> szintek; //vektor amiben minden szinthez minden ráta meg van adva
double mutrata;
double szamlalo; // ebbe fogom összeadni a rátákat
int szintszam;
int instanciak;
int instancia;
int gen;
int maxMut;//maximálisan begyűjthető mutációk száma
double randdriver;// mutrata*drivergének száma/összes   ~ 0.84*70/20000
double s_acd;
double s_scd;
double s_scdd; //fittnessz erősség
int sumnacd;  // összesített driver mutációk az adott szinten acd,scd,scdplusdből
int sumnscd;
int sumnscdplusd;
double utolsodelta;
double idosum;
double Nsuminstanc=0;
double Nnegyzetatlag_Ref=0;
double Nnegyzetatlag_Varied=0;
double supernegyzetart_Varied;
double Natlag_Ref=0;
double Natlag_Varied;
int vartosszes;
int Nrak;
double B;
int rakosodas;
double STOP;
double Dszum=0;
double Dszumnegyzet;
double Dnegyzetszum=0;
double Dnegyzetatlag=0;
double Dn_atlagsum=0;
double D_natlag=0;
int D_n=0;
int D_n_szamlalo=0;


ifstream f;

double superatl_Ref=0;
double superatl_Varied=0;
double summ;
double summNegyzet;
double sosszeg;

int Nmutmax=0;
char flag_optimal; //opt gamma?
char flag_stemMut; //őssejt mut?
char flag_gnuplot;
char flag_idKiir;
string outputDir;

int main()
{
    struct timeval time;
    gettimeofday(&time,NULL);

    double seed =(time.tv_sec * 100) + (time.tv_usec / 100);
    srand(seed);
    f.open("adatok.txt");
    /*cout.precision(17);
    cout<<endl<<seed<<endl<<endl;
    cout.precision(6);*/

    string tempRead;
    /*fancyRead(f,tempRead);
    instanciak=stoi(tempRead);
    fancyRead(f,tempRead);
    legyartando=stoi(tempRead);
    fancyRead(f,tempRead);
    rakosodas=stoi(tempRead);
    fancyRead(f,tempRead);
    STOP=stod(tempRead);
    fancyRead(f,tempRead);
    ido=stod(tempRead);
    fancyRead(f,tempRead);
    szintszam=stoi(tempRead);
    fancyRead(f,tempRead);
    mutrata=stod(tempRead);
    fancyRead(f,tempRead);
    s_acd=stod(tempRead);
    fancyRead(f,tempRead);
    s_scdd=stod(tempRead);
    fancyRead(f,tempRead);
    s_scd=stod(tempRead);
    fancyRead(f,tempRead);
    B=stod(tempRead);
    fancyRead(f,tempRead);
    maxMut=stoi(tempRead);
    fancyRead(f,tempRead);
    flag_optimal=tempRead.at(0);
    fancyRead(f,tempRead);
    flag_stemMut=tempRead.at(0);*/


    if(f.fail()){
        cout<<"error: adatok.txt missing"<<endl;
        return -1;
    }

    while(!f.eof()){

        f>>tempRead;
        if(tempRead=="instanciak="){
            f>>instanciak;
        }else if(tempRead=="legyartando="){
            f>>legyartando;
        }else if(tempRead=="ido="){
            f>>ido;
        }else if(tempRead=="rakosodas="){
            f>>rakosodas;
        }else if(tempRead=="STOP="){
            f>>STOP;
        }else if(tempRead=="szintszam="){
            f>>szintszam;
        }else if(tempRead=="mutrata="){
            f>>mutrata;
        }else if(tempRead=="s_acd="){
            f>>s_acd;
        }else if(tempRead=="s_scd="){
            f>>s_scd;
        }else if(tempRead=="s_scdd="){
            f>>s_scdd;
        }else if(tempRead=="B="){
            f>>B;
        }else if(tempRead=="maxMut="){
            f>>maxMut;
        }else if(tempRead=="flag_optimal="){
            f>>flag_optimal;
        }else if(tempRead=="flag_stemMut="){
            f>>flag_stemMut;
        }else if(tempRead=="outputDir="){
            f>>outputDir;
        }else if(tempRead=="flag_gnuplot="){
            f>>flag_gnuplot;
        }else if(tempRead=="flag_idKiir="){
            f>>flag_idKiir;
        }else if(tempRead=="matrix="){
            for(int i=0; i<szintszam+1;i++)
            {
                szint ujszint;
                f>>ujszint.delta;
                f>>ujszint.Nk;
                f>>ujszint.pk;
                //cout<<ujszint.delta<<"\t";
                //cout<<ujszint.Nk<<"\t";
                //cout<<ujszint.pk<<endl;
                ujszint.qk=0;
                ujszint.racd=0;
                ujszint.rscd=0;
                ujszint.rscdplusd=0;
                szintek.push_back(ujszint);
            }
        }else{
            getline(f,tempRead);
        }
    }
    if(outputDir=="default"){
        outputDir="";
    }
    f.close();

    ofstream g;
    //ofstream incidencia;
    //ofstream t("Nperidomatrix.txt");
    ofstream p(outputDir+"vegallapotok.txt");
    //ofstream tt("idoinstanciaatlagmatrix.txt");
    /*ofstream h(outputDir+"szintek.txt");
    ofstream hg(outputDir+"histog.txt");
    ofstream sz(outputDir+"szintadatok.txt");
    ofstream id;
*/
    /*
    cout<<instanciak<<endl;
    cout<<legyartando<<endl;
    cout<<ido<<endl;
    cout<<szintszam<<endl;
    cout<<mutrata<<endl;
    cout<<s_acd<<endl;
    cout<<s_scdd<<endl;
    cout<<s_scd<<endl;
    cout<<B<<endl;
    cout<<maxMut<<endl;
    cout<<flag_optimal<<endl;
    cout<<flag_stemMut<<endl;
    cout<<outputDir<<endl;
    */
    string swapString="inst="+to_string(instanciak);

    ifstream script("gnuplotscriptOriginal.plt");
    if(script.fail()){
        cout<<"error: gnuplotscriptOriginal.plt missing"<<endl;
        return -1;
    }
    txtsorcsere("inst=50",swapString,script,"gnuplotscript.plt");
    script.close();


    //double utolsodelta=1;

    vartosszes=0;

    for(int i=0;i<szintszam;i++)
    {
        vartosszes+=szintek[i].Nk;
    }


    script.open("gnuplotscript.plt");
    swapString="a="+to_string(vartosszes);
    txtsorcsere("a=301",swapString,script,"gnuplotscript.plt");
    script.close();



    szintek[0].pk=0;
    //ráták és qk kiszám
    szintek[0].rscd=0;
    szintek[0].rscdplusd=0;


   // summ=0;
   //summNegyzet=0;

    int Nsum;
    int sumnscd;
    int sumnacd;
    int sumnscdplusd;

    p<<setw(12)<<"inst";
    p<<setw(12)<<"Nsum";
    p<<setw(12)<<"racd";
    p<<setw(12)<<"rscd";
    p<<setw(12)<<"rscdd";
    p<<setw(12)<<"Vart";
    p<<setw(12)<<"Nsum/Vart";
    p<<setw(10)<<"D*";
    p<<setw(10)<<"D*szn"<<endl;

    supernegyzetart_Varied=0;
    superatl_Varied=0;
    Nnegyzetatlag_Varied=0;
    Natlag_Varied=0;
    Nrak=0;

for(instancia=0;instancia<instanciak;instancia++) // itt átlagolom ki egy forral a 100 adatból az átlagos osztódásokat az n. szinten (Dn függés y tengely)
{
    double atlagresz=0;
    double atlagresz_scd=0;
    double atlagresz_acd=0;
    double atlagresz_scdd=0;



sosszeg=(s_scd+s_scdd+s_acd);

    szintek[0].pk=0;



    if(flag_optimal=='y' || flag_optimal== 'Y') //Optimalis ratak
    {
         double gamma_optimalis;
         gamma_optimalis=pow(legyartando,-1.0/(double)szintszam);
         utolsodelta=1;


         if(pow(legyartando,-1.0/(double)szintszam)>0.5)
         {
            cout<<"Nem fog elmenni a megadott szintszamig ekkora N-el, mert az optimális gamma tul nagy lesz. Ha ez zavar, noveld meg N-t"<<endl;
         }

         for(int i=szintszam+1;i>0;i--)
         {

            if(i==szintszam)
            {
                szintek[i].delta=utolsodelta;


            }
            else if(i==szintszam-1){
                szintek[i].delta=1.0*gamma_optimalis;

            }
            else
            {
                szintek[i].delta=gamma_optimalis*szintek[i+1].delta;

            }
         }


         for(int i=0;i<szintszam+1;i++)
         {
            if(i==szintszam)
            {
                szintek[i].qk=1;
                szintek[i].pk=1;
                szintek[i].racd=(1.0/szintek[i].Nk)*2*szintek[i].delta*(1-szintek[i].pk)/(szintek[i].pk*szintek[i].qk);
                szintek[i].rscd=(1.0/szintek[i].Nk)*szintek[i].delta*(1-szintek[i].qk)/szintek[i].qk;
                szintek[i].rscdplusd=(1.0/szintek[i].Nk)*szintek[i].delta*(1/szintek[i].qk);
            }
            else
            {
                szintek[i].qk=2*(szintek[i].delta/szintek[i+1].delta)/szintek[i].pk;
                szintek[i].racd=(1.0/szintek[i].Nk)*szintek[i+1].delta*(1-szintek[i].pk);
                szintek[i].rscd=(1.0/szintek[i].Nk)*szintek[i].delta*(1-szintek[i].qk)/szintek[i].qk;
                szintek[i].rscdplusd=(1.0/szintek[i].Nk)*szintek[i].delta*(1/szintek[i].qk);

            }
        //szintek[i].racd=(1.0/szintek[i].Nk)*2*szintek[i].delta*(1-szintek[i].pk)/(szintek[i].pk*szintek[i].qk);

            if(szintek[i].pk>1 || szintek[i].qk>1)
            {
                cout<<"baj van, a pk vagy qk nagyobb lett 1-nel, valassz jobb pk-t ezen a szinten:"<<i<<endl;
                cout<<szintek[i].pk<<endl;
                cout<<szintek[i].qk<<endl;

                if(gamma_optimalis>0.5)
                {
                    cout<<"Az opt. gamma it túl nagy (>0.5). Próbáld meg nagyobb N-re ha több szintet akarsz (az optimális görbe csak idáig tart)";
                }

                return 1;
            }
                if((szintek[i].delta/szintek[i+1].delta)>0.5 && i!=szintszam)
                {
                    cout<<"baj van, a delta(k)/delta(k+1) nagyobb lett 0.5-nel ezen a szinten:"<<i<<endl;
                    return 3;
                }

         }



    }





    else if(flag_optimal=='n' || flag_optimal== 'N') // manualis ratak
    {
            for(int i=0; i<szintszam;i++)
            {
                szint ujszint;
                f>>ujszint.delta;
                f>>ujszint.Nk;
                f>>ujszint.pk;
                ujszint.qk=0;
                ujszint.racd=0;
                ujszint.rscd=0;
                ujszint.rscdplusd=0;
                szintek.push_back(ujszint);

            }
                //double utolsodelta=1;

            for(int i=0;i<szintszam;i++)
            {
                if(i==szintszam-1)
                {
                    szintek[i].qk=1;
                    szintek[i].pk=1;
                    szintek[i].racd=(1.0/szintek[i].Nk)*2*szintek[i].delta*(1-szintek[i].pk)/(szintek[i].pk*szintek[i].qk);
                    szintek[i].rscd=(1.0/szintek[i].Nk)*szintek[i].delta*(1-szintek[i].qk)/szintek[i].qk;
                    szintek[i].rscdplusd=(1.0/szintek[i].Nk)*szintek[i].delta*(1/szintek[i].qk);
                }
                else
                {
                    szintek[i].qk=2*(szintek[i].delta/szintek[i+1].delta)/szintek[i].pk;
                    szintek[i].racd=(1.0/szintek[i].Nk)*szintek[i+1].delta*(1-szintek[i].pk);
                    szintek[i].rscd=(1.0/szintek[i].Nk)*szintek[i].delta*(1-szintek[i].qk)/szintek[i].qk;
                    szintek[i].rscdplusd=(1.0/szintek[i].Nk)*szintek[i].delta*(1/szintek[i].qk);

                }
                //szintek[i].racd=(1.0/szintek[i].Nk)*2*szintek[i].delta*(1-szintek[i].pk)/(szintek[i].pk*szintek[i].qk);

                if(szintek[i].pk>1 || szintek[i].qk>1)
                {
                    cout<<"baj van, a pk vagy qk nagyobb lett 1-nel, valassz jobb pk-t ezen a szinten:"<<i<<endl;
                    cout<<szintek[i].pk<<endl;
                    cout<<szintek[i].qk<<endl;
                    return 1;
                }
                if((szintek[i].delta/szintek[i+1].delta)>0.5 && i!=szintszam-1)
                {
                    cout<<"baj van, a delta(k)/delta(k+1) nagyobb lett 0.5-nel ezen a szinten:"<<i<<endl;
                    return 3;
                }

            }



    }
     else
        {
            cout<<"Nem sikerult valasztani..."<<endl;
            return 1;
        }

    szintek[0].rscd=0;
    szintek[0].rscdplusd=0;



   /* string filenameRoot("Incidencia_");
    stringstream ss;
    ss << filenameRoot<<instancia<< ".txt";

    incidencia.open(outputDir+ss.str() );
*/

    string filenameRoot2("kimenoadatok_");
    stringstream sss;
    sss << filenameRoot2 <<instancia<< ".txt";

    g.open(outputDir+sss.str() );

  /*  string filenameRoot3("szintadatok_");
    stringstream ssss;
    ssss << filenameRoot3 <<instancia<< ".txt";

    string instString=to_string(instancia);
    sz.open(outputDir+ssss.str() );
    id.open(outputDir+"idKiir"+instString+".txt");
*/

    for(int i=0;i<szintszam+1;i++)
    {
        x[i].clear();

    }//inicializalja az adatszerkezetet



    //Poisson generálás a sima mutációknak

    double lambda=mutrata;
    vector<double> poi;
    int k=-1;
    int kfakt=1;

    do{
        k++;
        kfakt*=(k>0)? k : 1;
        double uj;
        uj=(1.0/(double)kfakt)*((double)pow(lambda,k))*exp(-lambda);
        //cout<<uj<<endl;
        poi.push_back(uj);
    }while(poi[k]>(1/(double)RAND_MAX));

    int idCounter=0;
    sejt ossejt;
    //ossejt properties beallitasa
    ossejt.oszt_szam=0;
    ossejt.mut_szam=0;
    ossejt.racd=szintek[0].racd;
    ossejt.rscd=szintek[0].rscd;
    ossejt.rscdplusd=szintek[0].rscdplusd;
    ossejt.n_acd=0;
    ossejt.n_scd=0;
    ossejt.n_scdplusd=0;
    ossejt.id=idCounter;


    x[0].push_back(ossejt);


    idoskala = 1.0/(ossejt.racd+ossejt.rscd+ossejt.rscdplusd);
    // Kezdetben csak egy sejt van úgyhogy annak a rátái határozzák meg az idõskálát a fenti módon, ha nem tévedek
    double eltelt_ido=0;
    double legyartott_sejtek=0;
    //double leszivas_szamlalo=0;
    /* g<<"#idő\t";
    for(int i=0;i<szintszam;i++)
    {
        g<<i<<".\t";
    }
    g<<endl;*/
    vector<double> szintRataSum;
    for(int i=0; i<szintszam;i++){
        szintRataSum.push_back(szintek[i].racd+szintek[i].rscd+szintek[i].rscdplusd);
    }

    gen=0;
    idosum=0;

    while((legyartott_sejtek < legyartando && Nsum<STOP*vartosszes) || x[szintszam-1].size()==0 )
    {


    for(int i=0; i<szintszam;i++)
    {
        for(unsigned int j=0;j<x[i].size();j++)
        {
            if(i==szintszam-1)
            {
                x[i][j].racd=szintek[i].racd;
                x[i][j].rscd=szintek[i].rscd*((1+B)*szintek[i].Nk)/(szintek[i].Nk+B*x[i].size());
                x[i][j].rscdplusd=szintek[i].rscdplusd;



            }

            else
            {
                x[i][j].racd=(szintek[i].racd+x[i][j].n_acd*s_acd*szintRataSum[i])*((1+B)*szintek[i+1].Nk)/(szintek[i+1].Nk+B*x[i+1].size());
                x[i][j].rscd=(szintek[i].rscd+x[i][j].n_scd*s_scd*szintRataSum[i])*((1+B)*szintek[i].Nk)/(szintek[i].Nk+B*x[i].size());
                x[i][j].rscdplusd=(szintek[i].rscdplusd+x[i][j].n_scdplusd*s_scdd*szintRataSum[i])*(szintek[i+1].Nk+B*szintek[i+1].Nk)/(szintek[i+1].Nk+B*x[i+1].size())*(szintek[i].Nk+B*x[i].size())/(szintek[i].Nk+B*szintek[i].Nk);

                //x[i][j].racd=(szintek[i].racd+x[i][j].n_acd*s_acd*szintek[i].racd)*((1+B)*szintek[i+1].Nk)/(szintek[i+1].Nk+B*x[i+1].size());
                //x[i][j].rscd=(szintek[i].rscd+x[i][j].n_scd*s_scd*szintek[i].rscd)*((1+B)*szintek[i].Nk)/(szintek[i].Nk+B*x[i].size());
                //x[i][j].rscdplusd=(szintek[i].rscdplusd+x[i][j].n_scdplusd*s_scdd*szintek[i].rscdplusd)*(szintek[i+1].Nk+B*szintek[i+1].Nk)/(szintek[i+1].Nk+B*x[i+1].size())*(szintek[i].Nk+B*x[i].size())/(szintek[i].Nk+B*szintek[i].Nk);

            }
                if(x[i][j].rscdplusd<0) x[i][j].rscdplusd=0;
                if(x[i][j].rscd<0) x[i][j].rscd=0;
                if(x[i][j].racd<0) x[i][j].racd=0;
        }


    }


    Nsum=0;
    sumnacd=0;
    sumnscd=0;
    sumnscdplusd=0;
/*
     if(gen==0)
     {
         incidencia<<setw(12)<<"#ido";
         incidencia<<setw(12)<<"#gen";
         incidencia<<setw(12)<<"#N";
         incidencia<<setw(12)<<"#nacd";
         incidencia<<setw(12)<<"#nscd";
         incidencia<<setw(12)<<"#nscdd";
         incidencia<<setw(12)<<"#Nidoatlag";
         incidencia<<setw(12)<<"#acdiatl";
         incidencia<<setw(12)<<"#scdiatl";
         incidencia<<setw(12)<<"#scddiatl";
         incidencia<<setw(12)<<"#Nmut";
         incidencia<<setw((maxMut+1)*12)<<"#Nmut_iatl";
         incidencia<<setw(12)<<endl;

         sz<<setw(12)<<"#N1";
         sz<<setw(12)<<"#N2";
         sz<<setw(12)<<"#N3";
         sz<<setw(12)<<"#N4";
         sz<<setw(12)<<"#n_scd1";
         sz<<setw(12)<<"#n_scd2";
         sz<<setw(12)<<"#n_scd3";
         sz<<setw(12)<<"#n_scd4"<<endl;


     }
*/
        // Kell minden körben egy randomszám ami eldönti, melyik esesmény fog bekövetkezni
        szamlalo=0;

        double melyik_esemeny;
        vector <lista> kommulativ_lista;
        for(int j=0;j<szintszam+1;j++) //ráták felhalmozása
        {

            for( unsigned int l=0;l<x[j].size();l++)
            {
                lista ujlistaelem;
                ujlistaelem.sejtszam=l;
                ujlistaelem.szintszam=j;

                szamlalo+=x[j][l].rscd;
                ujlistaelem.ertek=szamlalo;
                ujlistaelem.tipus=0;
                kommulativ_lista.push_back(ujlistaelem);

                szamlalo+=x[j][l].rscdplusd;
                ujlistaelem.tipus=1;
                ujlistaelem.ertek=szamlalo;
                kommulativ_lista.push_back(ujlistaelem);

                szamlalo+=x[j][l].racd;
                ujlistaelem.tipus=2;
                ujlistaelem.ertek=szamlalo;
                kommulativ_lista.push_back(ujlistaelem);
            }

        }


        double r=0;
        while(r==0)
        {
            r=(double)rand();
        }

        melyik_esemeny=(szamlalo*r/(double)RAND_MAX);

        int esemeny_szintje;
        int esemeny_sejtszama;
        int esemeny_tipus; // 0-rscd, 1 rscdplusd, 2 acd

        //Esemény helyének keresése felezéses módszerrel

        unsigned int kereso=((kommulativ_lista.size()/2));
        unsigned int keresomax=kommulativ_lista.size();
        unsigned int keresomin=0;
        //cout<<"\n\n\n kereso: "<<kereso<<"\n szamlalo: "<<szamlalo<<"\n kommlist size: "<<kommulativ_lista.size();
        //cout<<"\n melyik esem: "<<melyik_esemeny<<"n\n\n";
        if(melyik_esemeny==szamlalo)
        {
            kereso=kommulativ_lista.size();
        }
        else if(kereso!=1)
        {
            while(!(kommulativ_lista[kereso-1].ertek>=melyik_esemeny && kommulativ_lista[kereso-2].ertek<melyik_esemeny))
            {
                if(kommulativ_lista[kereso-1].ertek<melyik_esemeny)
                {
                    keresomin=kereso;
                    if((keresomax-keresomin)/2>0) {kereso=keresomin+(keresomax-keresomin)/2;}else {kereso++;}
                    if(kereso==kommulativ_lista.size()) {break;}
                }
                else
                {
                    keresomax=kereso;
                    if((keresomax-keresomin)/2>0) {kereso=keresomin+(keresomax-keresomin)/2;}else {kereso--;}
                    if(kereso==1) {break;}
                }
                //cout<<"max"<<keresomax<<endl;
                //cout<<keresomin<<endl;
                //cout<<kommulativ_lista[kereso-1].ertek<<endl<<"kell "<<melyik_esemeny<<endl;
            }
        }
        else //csak akkor lehet ilyenkor egy a kereső, ha csak az őssejt van a rendszerben (csúnya megoldás de már éjfél van)
        {
            kereso=3;
        }
        esemeny_sejtszama=kommulativ_lista[kereso-1].sejtszam;
        esemeny_szintje=kommulativ_lista[kereso-1].szintszam;
        esemeny_tipus=kommulativ_lista[kereso-1].tipus;

   /*     if(esemeny_szintje==0 && esemeny_tipus!=2)
        {
            cout<<"error: az őssejt rosszul osztódott"<<endl;
            cout<<kereso<<endl;
        }*/


        //felezéses módszer vége

        //itt bekövetkezik az esemény
        x[esemeny_szintje][esemeny_sejtszama].oszt_szam++;

        if(esemeny_tipus==0) //scd
        {
            idCounter++;
            sejt ujsejt;
            ujsejt.n_acd=x[esemeny_szintje][esemeny_sejtszama].n_acd;
            ujsejt.n_scd=x[esemeny_szintje][esemeny_sejtszama].n_scd;
            ujsejt.n_scdplusd=x[esemeny_szintje][esemeny_sejtszama].n_scdplusd;
            ujsejt.oszt_szam=x[esemeny_szintje][esemeny_sejtszama].oszt_szam;
            ujsejt.mut_szam=x[esemeny_szintje][esemeny_sejtszama].mut_szam;
            ujsejt.id=idCounter;

            double randmut=rand()/(double)RAND_MAX;
            double randpoi=0;


            for(int i=0;i<k;i++)
            {
                randpoi+=poi[i];
                if(randmut<=randpoi)
                {

                    for(int j=0; j<i;j++)
                    {
                        if(ujsejt.mut_szam<maxMut){
                            ujsejt.mut_szam++;
                            double randmut=rand()/(double)RAND_MAX;

                            if(randmut<=1.0/3.0)
                            {
                                ujsejt.n_acd++;
                                history newHistory;
                                newHistory.idopont=eltelt_ido;
                                newHistory.szint=esemeny_szintje;
                                newHistory.type="acd";
                                ujsejt.historyVector.push_back(newHistory);
                            }

                            if(1.0/3.0<randmut && randmut<2.0/3.0)
                            {
                                ujsejt.n_scd++;
                                history newHistory;
                                newHistory.idopont=eltelt_ido;
                                newHistory.szint=esemeny_szintje;
                                newHistory.type="scd";
                                ujsejt.historyVector.push_back(newHistory);
                            }

                            if(2.0/3.0<=randmut)
                            {
                                ujsejt.n_scdplusd++;
                                history newHistory;
                                newHistory.idopont=eltelt_ido;
                                newHistory.szint=esemeny_szintje;
                                newHistory.type="scdplusd";
                                ujsejt.historyVector.push_back(newHistory);

                            }


                        }

                    }


                    break;


                }
                 /*if(ujsejt.historyVector.size()>0){
                                cout<<endl<<"brekk pont: "<<ujsejt.historyVector.size()<<endl;
                             }*/

            }



            randmut=rand()/(double)RAND_MAX;
            randpoi=0;
            for(int i=0;i<k;i++)
            {
                randpoi+=poi[i];
                if(randmut<=randpoi)
                {

                    for(int j=0; j<i;j++)
                    {
                        if(x[esemeny_szintje][esemeny_sejtszama].mut_szam<maxMut){
                            x[esemeny_szintje][esemeny_sejtszama].mut_szam++;
                            double randmut=rand()/(double)RAND_MAX;

                            if(randmut<=1.0/3.0)
                            {
                                x[esemeny_szintje][esemeny_sejtszama].n_acd++;
                                history newHistory;
                                newHistory.idopont=eltelt_ido;
                                newHistory.szint=esemeny_szintje;
                                newHistory.type="acd";
                                x[esemeny_szintje][esemeny_sejtszama].historyVector.push_back(newHistory);
                            }

                            if(1.0/3.0<randmut && randmut<2.0/3.0)
                            {
                                x[esemeny_szintje][esemeny_sejtszama].n_scd++;
                                history newHistory;
                                newHistory.idopont=eltelt_ido;
                                newHistory.szint=esemeny_szintje;
                                newHistory.type="scd";
                                x[esemeny_szintje][esemeny_sejtszama].historyVector.push_back(newHistory);
                            }

                            if(2.0/3.0<=randmut)
                            {
                                x[esemeny_szintje][esemeny_sejtszama].n_scdplusd++;
                                history newHistory;
                                newHistory.idopont=eltelt_ido;
                                newHistory.szint=esemeny_szintje;
                                newHistory.type="scdplusd";
                                x[esemeny_szintje][esemeny_sejtszama].historyVector.push_back(newHistory);

                            }
                        }
                    }
                    break;
                }


            }


            x[esemeny_szintje].push_back(ujsejt);


        }
        else if(esemeny_tipus==1) //scdplusd
        {
            if(esemeny_szintje!=szintszam-1)
            {
                sejt ujsejt;
                ujsejt.oszt_szam=x[esemeny_szintje][esemeny_sejtszama].oszt_szam;
                ujsejt.mut_szam=x[esemeny_szintje][esemeny_sejtszama].mut_szam;
                ujsejt.n_acd=x[esemeny_szintje][esemeny_sejtszama].n_acd;
                ujsejt.n_scd=x[esemeny_szintje][esemeny_sejtszama].n_scd;
                ujsejt.n_scdplusd=x[esemeny_szintje][esemeny_sejtszama].n_scdplusd;
                ujsejt.id=x[esemeny_szintje][esemeny_sejtszama].id;

                idCounter++;
                sejt ujsejt2;
                ujsejt2.oszt_szam=x[esemeny_szintje][esemeny_sejtszama].oszt_szam;
                ujsejt2.mut_szam=x[esemeny_szintje][esemeny_sejtszama].mut_szam;
                ujsejt2.n_acd=x[esemeny_szintje][esemeny_sejtszama].n_acd;
                ujsejt2.n_scd=x[esemeny_szintje][esemeny_sejtszama].n_scd;
                ujsejt2.n_scdplusd=x[esemeny_szintje][esemeny_sejtszama].n_scdplusd;
                ujsejt2.id=idCounter;

                double randmut=rand()/(double)RAND_MAX;
                //cout<<randmut<<endl;
                double randpoi=0;
                for(int i=0;i<k;i++)
                {
                    randpoi+=poi[i];
                    if(randmut<=randpoi)
                    {


                        for(int j=0; j<i;j++)
                        {
                            if(ujsejt.mut_szam<maxMut){
                                ujsejt.mut_szam++;
                                double randmut=rand()/(double)RAND_MAX;

                                if(randmut<=1.0/3.0)
                                {
                                    ujsejt.n_acd++;
                                    history newHistory;
                                    newHistory.idopont=eltelt_ido;
                                    newHistory.szint=esemeny_szintje;
                                    newHistory.type="acd";
                                    ujsejt.historyVector.push_back(newHistory);
                                }

                                if(1.0/3.0<randmut && randmut<2.0/3.0)
                                {
                                    ujsejt.n_scd++;
                                    history newHistory;
                                    newHistory.idopont=eltelt_ido;
                                    newHistory.szint=esemeny_szintje;
                                    newHistory.type="scd";
                                    ujsejt.historyVector.push_back(newHistory);
                                }

                                if(2.0/3.0<=randmut)
                                {
                                    ujsejt.n_scdplusd++;
                                    history newHistory;
                                    newHistory.idopont=eltelt_ido;
                                    newHistory.szint=esemeny_szintje;
                                    newHistory.type="scdplusd";
                                    ujsejt.historyVector.push_back(newHistory);
                                }
                            }
                        }
                        break;
                    }
                }




                randmut=rand()/(double)RAND_MAX;
                randpoi=0;
               for(int i=0;i<k;i++)
                {
                    randpoi+=poi[i];
                    if(randmut<=randpoi)
                    {


                        for(int j=0; j<i;j++)
                        {
                            if(ujsejt2.mut_szam<maxMut){
                                ujsejt2.mut_szam++;
                                double randmut=rand()/(double)RAND_MAX;

                                if(randmut<=1.0/3.0)
                                {
                                    ujsejt2.n_acd++;
                                    history newHistory;
                                    newHistory.idopont=eltelt_ido;
                                    newHistory.szint=esemeny_szintje;
                                    newHistory.type="acd";
                                    ujsejt2.historyVector.push_back(newHistory);
                                }

                                if(1.0/3.0<randmut && randmut<2.0/3.0)
                                {
                                    ujsejt2.n_scd++;
                                    history newHistory;
                                    newHistory.idopont=eltelt_ido;
                                    newHistory.szint=esemeny_szintje;
                                    newHistory.type="scd";
                                    ujsejt2.historyVector.push_back(newHistory);
                                }

                                if(2.0/3.0<=randmut)
                                {
                                    ujsejt2.n_scdplusd++;
                                    history newHistory;
                                    newHistory.idopont=eltelt_ido;
                                    newHistory.szint=esemeny_szintje;
                                    newHistory.type="scdplusd";
                                    ujsejt2.historyVector.push_back(newHistory);
                                }
                            }

                          }
                          break;
                    }
                }




                x[esemeny_szintje+1].push_back(ujsejt);
                x[esemeny_szintje+1].push_back(ujsejt2);


            }
            else
            {
                legyartott_sejtek+=2;
                if(legyartott_sejtek == legyartando-2 || legyartott_sejtek == legyartando-1 ){
                    D_n += x[esemeny_szintje][esemeny_sejtszama].oszt_szam;
                    D_n_szamlalo++;

                }
            }


            std::vector<sejt>::iterator it=x[esemeny_szintje].begin();
            for(int i=0;i<esemeny_sejtszama;i++)it++;
            x[esemeny_szintje].erase(it);
        }
        else //acd
        {
            if(esemeny_szintje!=(szintszam-1))
            {
                idCounter++;
                sejt ujsejt;
                ujsejt.oszt_szam=x[esemeny_szintje][esemeny_sejtszama].oszt_szam;
                ujsejt.mut_szam=x[esemeny_szintje][esemeny_sejtszama].mut_szam;
                ujsejt.n_acd=x[esemeny_szintje][esemeny_sejtszama].n_acd;
                ujsejt.n_scd=x[esemeny_szintje][esemeny_sejtszama].n_scd;
                ujsejt.n_scdplusd=x[esemeny_szintje][esemeny_sejtszama].n_scdplusd;
                ujsejt.id=idCounter;

                double randmut=rand()/(double)RAND_MAX;
                double randpoi=0;
                for(int i=0;i<k;i++)
                {
                    randpoi+=poi[i];
                    if(randmut<=randpoi)
                    {


                    for(int j=0; j<i;j++)
                    {
                        if(ujsejt.mut_szam<maxMut){
                            ujsejt.mut_szam++;
                            double randmut=rand()/(double)RAND_MAX;

                            if(randmut<=1.0/3.0)
                            {
                                ujsejt.n_acd++;
                                history newHistory;
                                newHistory.idopont=eltelt_ido;
                                newHistory.szint=esemeny_szintje;
                                newHistory.type="acd";
                                ujsejt.historyVector.push_back(newHistory);
                            }

                            if(1.0/3.0<randmut && randmut<2.0/3.0)
                            {
                                ujsejt.n_scd++;
                                history newHistory;
                                newHistory.idopont=eltelt_ido;
                                newHistory.szint=esemeny_szintje;
                                newHistory.type="scd";
                                ujsejt.historyVector.push_back(newHistory);
                            }

                            if(2.0/3.0<=randmut)
                            {
                                ujsejt.n_scdplusd++;
                                history newHistory;
                                newHistory.idopont=eltelt_ido;
                                newHistory.szint=esemeny_szintje;
                                newHistory.type="scdplusd";
                                ujsejt.historyVector.push_back(newHistory);
                            }
                        }
                    }
                    break;
                    }
                }

                randmut=rand()/(double)RAND_MAX;
                randpoi=0;


              if(esemeny_szintje==0)
              {

                if(flag_stemMut =='y' || flag_stemMut == 'Y')
                {
                      for(int i=0;i<k;i++)
                      {
                        randpoi+=poi[i];
                        if(randmut<=randpoi)
                        {


                            for(int j=0; j<i;j++)
                            {
                                if(x[esemeny_szintje][esemeny_sejtszama].mut_szam<maxMut){
                                    x[esemeny_szintje][esemeny_sejtszama].mut_szam++;
                                    double randmut=rand()/(double)RAND_MAX;

                                    if(randmut<=1.0/3.0)
                                    {
                                        x[esemeny_szintje][esemeny_sejtszama].n_acd++;
                                        history newHistory;
                                        newHistory.idopont=eltelt_ido;
                                        newHistory.szint=esemeny_szintje;
                                        newHistory.type="acd";
                                        ujsejt.historyVector.push_back(newHistory);

                                    }

                                    if(1.0/3.0<randmut && randmut<2.0/3.0)
                                    {
                                        x[esemeny_szintje][esemeny_sejtszama].n_scd++;
                                        history newHistory;
                                        newHistory.idopont=eltelt_ido;
                                        newHistory.szint=esemeny_szintje;
                                        newHistory.type="scd";
                                        ujsejt.historyVector.push_back(newHistory);
                                    }

                                    if(2.0/3.0<=randmut)
                                    {
                                        x[esemeny_szintje][esemeny_sejtszama].n_scdplusd++;
                                        history newHistory;
                                        newHistory.idopont=eltelt_ido;
                                        newHistory.szint=esemeny_szintje;
                                        newHistory.type="scdplusd";
                                        ujsejt.historyVector.push_back(newHistory);
                                    }
                                }
                            }
                            break;
                        }
                      }

                }
              }

              else
              {
                 for(int i=0;i<k;i++)
                      {
                        randpoi+=poi[i];
                        if(randmut<=randpoi)
                        {


                            for(int j=0; j<i;j++)
                            {
                                if(x[esemeny_szintje][esemeny_sejtszama].mut_szam<maxMut){
                                    x[esemeny_szintje][esemeny_sejtszama].mut_szam++;
                                    double randmut=rand()/(double)RAND_MAX;

                                    if(randmut<=1.0/3.0)
                                    {
                                        x[esemeny_szintje][esemeny_sejtszama].n_acd++;
                                        history newHistory;
                                        newHistory.idopont=eltelt_ido;
                                        newHistory.szint=esemeny_szintje;
                                        newHistory.type="acd";
                                        x[esemeny_szintje][esemeny_sejtszama].historyVector.push_back(newHistory);
                                    }

                                    if(1.0/3.0<randmut && randmut<2.0/3.0)
                                    {
                                        x[esemeny_szintje][esemeny_sejtszama].n_scd++;
                                        history newHistory;
                                        newHistory.idopont=eltelt_ido;
                                        newHistory.szint=esemeny_szintje;
                                        newHistory.type="scd";
                                        x[esemeny_szintje][esemeny_sejtszama].historyVector.push_back(newHistory);
                                    }

                                    if(2.0/3.0<=randmut)
                                    {
                                        x[esemeny_szintje][esemeny_sejtszama].n_scdplusd++;
                                        history newHistory;
                                        newHistory.idopont=eltelt_ido;
                                        newHistory.szint=esemeny_szintje;
                                        newHistory.type="scdplusd";
                                        x[esemeny_szintje][esemeny_sejtszama].historyVector.push_back(newHistory);
                                    }
                                }
                            }
                            break;
                        }
                      }
              }





                x[esemeny_szintje+1].push_back(ujsejt);


            }
            else
            {
                legyartott_sejtek++;
                if(legyartott_sejtek==legyartando-1){
                    D_n += x[esemeny_szintje][esemeny_sejtszama].oszt_szam;
                    D_n_szamlalo++;
                }

            }







        }

            vector<int> Nmut;
            for(int i=0;i<maxMut+1;i++){
                Nmut.push_back(0);
            }
            for(unsigned int i=0;i<x.size();i++)
            {
                Nsum+=x[i].size();

                for(unsigned int j=0; j<x[i].size();j++)
                {
                    sumnacd+=x[i][j].n_acd;
                    sumnscd+=x[i][j].n_scd;
                    sumnscdplusd+=x[i][j].n_scdplusd;

                    if(x[i][j].mut_szam<Nmut.size()-1){
                        Nmut[x[i][j].mut_szam]++;
                    }
                    else{
                        for(unsigned int k=0;k<(x[i][j].mut_szam-(Nmut.size()-1));k++){
                            Nmut.push_back(0);
                        }
                        Nmut[x[i][j].mut_szam]++;
                    }
                    if(x[i][j].mut_szam>Nmutmax){
                        Nmutmax=x[i][j].mut_szam;
                    }
                }
            }

/*
            int szintsumnscd1=0;
            int szintsumnscd2=0;
            int szintsumnscd3=0;
            int szintsumnscd4=0;

            for(unsigned int j=0; j<x[1].size();j++)
                {
                    szintsumnscd1+=x[1][j].n_scd;
                }

            for(unsigned int j=0; j<x[2].size();j++)
                {
                    szintsumnscd2+=x[2][j].n_scd;
                }

            for(unsigned int j=0; j<x[3].size();j++)
                {
                    szintsumnscd3+=x[3][j].n_scd;
                }

            for(unsigned int j=0; j<x[4].size();j++)
                {
                    szintsumnscd4+=x[4][j].n_scd;
                }
*/
            vector<double> Nmut_atlagresz;
            for(int i=0;i<50;i++)
            {
                Nmut_atlagresz.push_back(0);
            }
            atlagresz+=(Nsum*idoskala);
            atlagresz_acd+=(sumnacd*idoskala);
            atlagresz_scd+=(sumnscd*idoskala);
            atlagresz_scdd+=(sumnscdplusd*idoskala);

            for(unsigned int i=0;i<Nmut.size();i++)
            {
                Nmut_atlagresz[i]+=Nmut[i]*idoskala;
            }

	/*
            incidencia<<setw(12)<<eltelt_ido;
            incidencia<<setw(12)<<gen;
            incidencia<<setw(12)<<Nsum;
            incidencia<<setw(12)<<sumnacd;
            incidencia<<setw(12)<<sumnscd;
            incidencia<<setw(12)<<sumnscdplusd;
            incidencia<<setw(12)<<(atlagresz/eltelt_ido);
            incidencia<<setw(12)<<(atlagresz_acd/eltelt_ido);
            incidencia<<setw(12)<<(atlagresz_scd/eltelt_ido);
            incidencia<<setw(12)<<(atlagresz_scdd/eltelt_ido);
            for(unsigned int i=0;i<maxMut+1;i++){
                incidencia<<setw(12)<<Nmut[i];
            }
            for(unsigned int i=0;i<maxMut+1;i++)
            {
                incidencia<<setw(12)<<Nmut_atlagresz[i]/eltelt_ido;
            }
            incidencia<<setw(12)<<endl;

            sz<<setw(12)<<x[1].size();
            sz<<setw(12)<<x[2].size();
            sz<<setw(12)<<x[3].size();
            sz<<setw(12)<<x[4].size();
            sz<<setw(12)<<szintsumnscd1;
            sz<<setw(12)<<szintsumnscd2;
            sz<<setw(12)<<szintsumnscd3;
            sz<<setw(12)<<szintsumnscd4<<endl;


*/

            idoskala=1.0/szamlalo;
            eltelt_ido+=idoskala;


            gen++;



} //eddig tart legyartani az osszes sejtet


    if(Nsum>vartosszes*rakosodas)
    {
        Nrak++;
    }


     Dn_atlagsum=0;
     Dnegyzetszum=0;
     p<<setw(12)<<instancia;
     p<<setw(12)<<Nsum;
     p<<setw(12)<<sumnacd;
     p<<setw(12)<<sumnscd;
     p<<setw(12)<<sumnscdplusd;
     p<<setw(12)<<vartosszes;
     p<<setw(12)<<100*Nsum/vartosszes<<"%";
     for(unsigned int j=0; j<x[szintszam-1].size();j++)
     {
          Dn_atlagsum+=x[szintszam-1][j].oszt_szam;
          Dnegyzetszum+=pow(x[szintszam-1][j].oszt_szam,2);
     }
     D_natlag=Dn_atlagsum/x[szintszam-1].size();
     Dszum+=D_natlag;
     Dnegyzetatlag=Dnegyzetszum/x[szintszam-1].size();
     Dszumnegyzet+=Dnegyzetatlag;
     p<<setw(10)<<D_natlag;
     p<<setw(10)<<Dnegyzetatlag<<endl;

        Nnegyzetatlag_Varied+=pow(Nsum,2);
        Natlag_Varied+=Nsum;





    //g.open("kimenoadatok.txt");
    //szint adatok kiiratas
    int w=17;
    g<<"Az egyes szintek adatai:"<<endl<<endl;
    g<<setw(5)<<"szint";
    g<<setw(8)<<"Nk";
    g<<setw(w)<<"racd";
    g<<setw(w)<<"rscd";
    g<<setw(w)<<"rscdplusd";
    g<<setw(w)<<"pk";
    g<<setw(w)<<"qk";
    g<<setw(w)<<"delta"<<endl;
    for(int i=0;i<szintszam+1;i++)
    {
        g<<setw(5)<<i;
        g<<setw(8)<<szintek[i].Nk;
        g<<setw(w)<<szintek[i].racd;
        g<<setw(w)<<szintek[i].rscd;
        g<<setw(w)<<szintek[i].rscdplusd;
        g<<setw(w)<<szintek[i].pk;
        g<<setw(w)<<szintek[i].qk;
        g<<setw(w)<<szintek[i].delta<<endl;
    }
    g<<endl;
    g<<setw(w)<<"s_scd, s_acd , s_scdd: "<<s_scd<<","<<s_acd<<","<<s_scdd<<endl;
    g.precision(17);
    g<<"Random seed:"<<seed<<endl;
    g.precision(6);
    g<<endl<<"A szimuláció alatt legyártott sejtek száma: "<<legyartott_sejtek<<endl;
    g<<endl<<"A szimuláció alatt eltelt idő: "<<eltelt_ido<<endl;


    for(int i=0; i<szintszam+1;i++)
    {
        g<<endl<<i<<". szint"<<endl;
        int sum=0;
        g<<"sejtek osztódási számai:"<<endl;
        for(unsigned int j=0; j<x[i].size();j++)
        {
                g<<x[i][j].oszt_szam<<" ";
                sum+=x[i][j].oszt_szam;
        }
        double atl=(double)sum/(double)(x[i].size());
        g<<endl<<"az átlagos osztódási szám: "<<atl<<endl;
        g<<"sejtek száma a szinten: ";
        g<<x[i].size()<<endl;
        g<<"ez "<<100*((double)x[i].size()/(double)szintek[i].Nk)<<"%-a a vártnak (Nk)"<<endl;

        sum=0;
        g<<"sejtek mutációs számai:"<<endl;
        for(unsigned int j=0; j<x[i].size();j++)
        {
                g<<x[i][j].mut_szam<<" ";
                sum+=x[i][j].mut_szam;


        }

        atl=(double)sum/(double)(x[i].size());
        g<<endl<<"az átlagos mutációs szám: "<<atl<<endl;

        sumnacd=0;
        g<<" acd driverek száma:"<<endl;
        for(unsigned int j=0; j<x[i].size();j++)
        {
                g<<x[i][j].n_acd<<" ";
                sumnacd+=x[i][j].n_acd;

        }

        double atlnacd=(double)sumnacd/(double)(x[i].size());
        g<<endl<<"az átlagos n_acd szám: "<<atlnacd<<endl;

        sumnscd=0;
        g<<" scd driverek száma:"<<endl;
        for(unsigned int j=0; j<x[i].size();j++)
        {
                g<<x[i][j].n_scd<<" ";
                sumnscd+=x[i][j].n_scd;

        }

        double atlnscd=(double)sumnscd/(double)(x[i].size());
        g<<endl<<"az átlagos n_scd szám: "<<atlnscd<<endl;

        sumnscdplusd=0;

        g<<" scdplusd driverek száma:"<<endl;
        for(unsigned int j=0; j<x[i].size();j++)
        {
                g<<x[i][j].n_scdplusd<<" ";
                sumnscdplusd+=x[i][j].n_scdplusd;

        }

        double atlnscdplusd=(double)sumnscdplusd/(double)(x[i].size());
        g<<endl<<"az átlagos n_scdplusd szám: "<<atlnscdplusd<<endl;
    }

if(flag_idKiir=='y'){
    /*id<<setw(15)<<"id";
    id<<setw(15)<<"tipus";
    id<<setw(15)<<"szint";
    id<<setw(15)<<"idopont";
    id<<endl;*/
    for(int i=0;i<szintszam;i++){
        for(int j=0;j<x[i].size();j++){
            for(int k=0;k<x[i][j].historyVector.size();k++){
               /* id<<setw(15)<<x[i][j].id;
                id<<setw(15)<<x[i][j].historyVector[k].type;
                id<<setw(15)<<x[i][j].historyVector[k].szint;
                id<<setw(15)<<x[i][j].historyVector[k].idopont;
                id<<endl;*/

            }

        }
        //id<<i<<endl;
    }
}

//ss.str().clear();
sss.str().clear();
/*ssss.str().clear();
incidencia.close();
*/g.close();
/*sz.close();
id.close();
*/
}//eddig tart az instancia

superatl_Varied=Natlag_Varied/instanciak;
supernegyzetart_Varied=Nnegyzetatlag_Varied/instanciak;

double Dsupernegyzetatlag=0;
double Dsuperatlag=0;

Dsupernegyzetatlag=Dszumnegyzet/instanciak;
Dsuperatlag=Dszum/instanciak;

p<<setw(8)<<"Ennyi sejt lett az összes szimuláció alatt a rendszerben a várthoz képest (Variált): "<<100*superatl_Varied/vartosszes<<"%"<<endl;
p<<setw(8)<<"Ennyi a sejtek számának szórása(Variált):"<< sqrt(supernegyzetart_Varied-pow(superatl_Varied,2))<<endl;
p<<setw(8)<<"A szórás a rendszer várt méretéhez képest(Variált): " << 100* (sqrt(supernegyzetart_Varied-pow(superatl_Varied,2)))/vartosszes<<"%"<<endl;
p<<setw(8)<<"n= "<<szintszam<<endl;
p<<setw(8)<<"D= "<< Dsuperatlag<<endl;
p<<setw(8)<<"D_n= "<< D_n/(double)instanciak<<endl;
p<<setw(8)<<"D_n_szamlalo= "<< D_n_szamlalo/(double)instanciak<<endl;
p<<setw(8)<<"szoras= "<<sqrt((Dsupernegyzetatlag)-pow((Dsuperatlag),2))<<endl;
p<<setw(8)<<"Nrak: "<<Nrak<<endl;

p.close();


script.open("gnuplotscript.plt");
swapString="maxmutszam="+to_string(Nmutmax  );
txtsorcsere("maxmutszam=0",swapString,script,"gnuplotscript.plt");
script.close();
if(flag_gnuplot=='y'){
    cout<<endl<<"Ábrák készítése..."<<endl;
    system("gnuplot gnuplotscript.plt");
    //system("gnuplot idoatlscript.plt");
    //system("gnuplot konvergenciaido.plt");
    cout<<"Ábrák elkészítve"<<endl;
}
    return 0;
}

void txtsorcsere(string ezt, string erre, ifstream& in,const char *inStreamFileName)
{
    ofstream out("tmp.txt");
    string tmp;
    while(true){
        getline(in,tmp);
        if(tmp==ezt){
            out<<erre<<endl;
        }else{
            out<<tmp<<endl;
        }
        if(in.eof()) break;
    }
    rename("tmp.txt",inStreamFileName);
    out.close();


}
void fancyRead(ifstream & ifs, string & tmp){

    while(true){
        getline(ifs,tmp);
        if(ifs.eof()){
            break;
        }else if(tmp.length()!=0){
            if(tmp.at(0)!='#') break;
        }
    }
}
