// g++ main_nographics.cpp colonies_univ.cpp dynamics_rods_3d.cpp -w -lpsapi -O3 -o colonies.exe

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
using namespace std;
#define SQR(x)  (x)*(x)
#include "functions.h"
#include "classes.h"

extern vector <Bacterium*> bact ;
extern double total_y_shift ;
extern ofstream runs ;
int sample ;

void err(char *reason)
{
  cout <<reason<<endl ; 
#ifndef __linux
  system("pause") ;
#endif  
  exit(0) ;
}

void err(char *reason, int a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifndef __linux
  system("pause") ;
#endif    
  exit(0) ;
}

void err(char *reason, double a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifndef __linux
  system("pause") ;
#endif    
  exit(0) ;
}

extern long long int xdr48 ;
extern long long tinit ; // real time at start

int main(int argc, char *argv[])
{
  NUM=new char[256] ; strcpy(NUM,"test_dir") ;
  int nsam=1,stopcond=0 ;
  start_file=NULL ;
  float max_y_shift,max_time ;
  int max_size=0, stc ;
  if (argc<3) { err(" Error:: at least two arguments needed, option and value, for example:\n"
                    "n [directory name]\n"
                    "s [fitness advantage]\n" 
                    "r [RNG seed]\n"
                    "N [no. of samples]\n"
                    "S [stop when: 1=only one type in the population]\n"
                    "\t[2 = after travelling [distance] (next argument)\n"
                    "\t[4 = after reaching size (next argument)\n"
										"\t[8 = after given time (next arg)]. Conditions can be combined e.g. S 1 S 8 1000\n"                    
                    "i [start filename]\n"
                    "e [eat rate]\n"
                    "E [elastic modulus]\n"
                    "w or h [width or height]\n"
                    "D [diffusion rate]\n"
                    "c0 [init. nutrient concentration]\n"
                    "z [zeta friction coefficient]\n"
                    "P [period for the sine at the bottom]\n"
                    "A [amplitude of the size at the botoom]\n"
                    "Program terminated. \n"); } 
  else { 
    for (int i=1;i<argc;i++) {
      switch (argv[i][0]) {
        case 'c': c0=atof(argv[++i]) ; cout <<"c0="<<c0<<" " ; break ;
        case 'D': diff_rate=atof(argv[++i]) ; cout <<"diff_rate="<<diff_rate<<" " ; break ;
        case 'e': Rod_shaped::eat_rate=atof(argv[++i]) ; cout <<"eat_rate="<<Rod_shaped::eat_rate<<" " ; break ;
        case 'E': Rod_shaped::E_bact=atof(argv[++i]) ; cout <<"E_bact="<<Rod_shaped::E_bact<<" " ; break ;
        case 'i': start_file=argv[++i] ; cout <<"startfile="<<start_file<<" " ;break ;
        case 'h': height=atoi(argv[++i]) ; cout <<"height="<<height<<" " ; break ;
        case 'n': NUM=argv[++i] ; cout <<"NUM="<<NUM<<" " ; break ;
        case 'N': nsam=atoi(argv[++i]) ; cout <<"nsam="<<nsam<<" " ; break ;
        case 'r': RAND=atoi(argv[++i]) ; cout <<"RAND="<<RAND<<" " ; break ;
        case 's': fit_adv=atof(argv[++i]) ; cout <<"s="<<fit_adv<<" " ; break ;
        case 'S': stc=atoi(argv[++i]) ; cout <<"stopcond="<<stc<<" " ; 
          if (stc==2) max_y_shift=atof(argv[++i]) ;
          if (stc==4) max_size=atoi(argv[++i]) ;
          if (stc==8) max_time=atof(argv[++i]) ;
          stopcond+=stc ; cout<<"total stopcond="<<stopcond<<" " ;
          break ;
        case 'w': width=atoi(argv[++i]) ; cout <<"width="<<width<<" " ; break ;
        case 'z': Rod_shaped::zeta=atof(argv[++i]) ; cout <<"zeta="<<Rod_shaped::zeta<<" " ; break ;
        case 'P': bottom_period=atof(argv[++i]) ; cout <<"bottom_period="<<bottom_period<<" " ; break ;
        case 'A': bottom_ampl=atof(argv[++i]) ; cout <<"bottom_ampl="<<bottom_ampl<<" " ; break ;
        default: err("unrecognized option.") ; break ;
      }
    }
    cout <<endl ;
//    NUM=argv[1] ;
//    nsam=atoi(argv[2]) ;
//    RAND=atoi(argv[3]) ;
  }
//  cout <<NUM<<" "<<" "<<nsam<<" "<<RAND<<endl ;
  _srand48(RAND) ;
  init();
  
  for (sample=0;sample<nsam;sample++) {
    if (sample>0) { // reset if necessary for more than 1 sample
      reset(false) ; 
    } 
    for (long long int i=0;;i++) {
      run() ;
      if (i%1000==0) {
        cout <<1.*(clock()-tinit)/CLOCKS_PER_SEC<<"\ts"<<sample<<": N="<<bact.size()<<" T="<<tt<<" #mut="<<count_types(1)<<endl ; 
        print_boxes() ;
      }
      if ((stopcond&1) && mut_occurred && is_only_one_type()==true) {
        char txt[256] ;
        sprintf(txt,"%s/runs_%d.dat",NUM,RAND) ; runs.open(txt,ios::app) ; 
        runs << "fixed: "<<bact[0]->type<<" " ; runs.close() ;  
        break ;
      }
      if ((stopcond&2) && -total_y_shift>=max_y_shift) {
        char txt[256] ;
        sprintf(txt,"%s/runs_%d.dat",NUM,RAND) ; runs.open(txt,ios::app) ; 
        runs << "max_y_reached" ; runs.close() ;  
        break ;
      }
      if ((stopcond&4) && bact.size()>=max_size) {
        char txt[256] ;
        sprintf(txt,"%s/runs_%d.dat",NUM,RAND) ; runs.open(txt,ios::app) ; 
        runs << "max_size_reached" ; runs.close() ;  
        break ;
      }
      if ((stopcond&8) && tt>=max_time) {
        char txt[256] ;
        sprintf(txt,"%s/runs_%d.dat",NUM,RAND) ; runs.open(txt,ios::app) ; 
      	runs <<"max_time reached"; runs.close() ;
        break ;
      }
    }
  }
	return 0 ;
}

void Rod_shaped::draw()
{    
}


