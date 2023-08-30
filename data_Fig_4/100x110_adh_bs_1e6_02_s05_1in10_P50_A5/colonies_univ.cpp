#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
using namespace std;
#include "functions.h"
#include "vecs.h"
#include "classes.h"
#include <time.h>

#define __MAIN
#include "params.h"


// -----------variables used in the main loop--------------------
// cc[][] is the concentration of food, dcc[][] is used to increment cc[][] in
// the procedure "diffusion"
// tt is time

char *start_file ;
double dt ;
double d0=2*r0 ; // width of the cell
double dx0 ; // diffusion grid length
char *NUM ;
int RAND=1 ;
int max_size=10000000 ;
int gr_ok=1 ; // if 1 then growth is simulated
bool mut_occurred ;
double tt ;
double **cc=NULL,**dcc=NULL,**sinkcc,**slime=NULL,*ccone=NULL,*dccone=NULL ;
#ifdef RECTBOX
  int **bulkbact=NULL ; // bulkbact[y][x] contains numbers of bacteria that are beyond the growing layer in a (x,y) box 1um x 1um
#endif
int tabx,taby, heightone;
float *front_thickness, min_front_thickness, mean_front_thickness, mean_front_roughness ;

double total_y_shift, prev_front_pos ;
double ycm ; // center of mass in y-direction
double yfront ; // position of the movign front

int ttt;  // number of simulation steps so far
int ttt_last_save, ttt_last_dt_update ; 
int how_many_steps_save ; 
int no_of_mutations=0 ;
bool reassign=false ;


vector <Bacterium*> bact ; // vector of cells

// these we don't save:

long long tinit ; // real time at start
int noncriterr ; // if set, it means that a non-critical error has occured and the last saved state has to be loaded 
double zmin,zmax,xmin,xmax,ymin,ymax ;  
int nointeracting ;

#ifdef CLOSEST
int noclosest ;
double avclosest ; // average number of closest neightbours
int nbox, boxwidth=int(10*d0) ;
vector <vector <vector <int> > > boxes ; // table of boxes (vectors) holding the cells
#else
int nbox=128, boxwidth=int(10*d0) ;
vector <vector <vector <int> > > boxes ;  // table of boxes (vectors) holding the cells
#endif

ofstream secs, params, runs ;
//#ifdef JUNCTIONS
vector <Junction*> junctions ;
//#endif


#ifdef __linux
#include <unistd.h>
typedef unsigned int DWORD ;
int memory_taken() // return memory available in MB
{
  long long int rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return (size_t)0L;		/* Can't open? */
	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
	{
		fclose( fp );
		return (size_t)0L;		/* Can't read? */
	}
	fclose( fp );
	long long int mem=((size_t)rss * (size_t)sysconf( _SC_PAGESIZE)) ;
	return (int) (mem/(1<<20));
}
#include <sys/sysinfo.h>
unsigned int freemem() // returns available memory in MB
{
  struct sysinfo s ;
  sysinfo(&s) ;
  return ((s.freeram)>>20) ;
}
#else
#include <windows.h>
#include <psapi.h>
int memory_taken() // return memory available in MB
{
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (int) (info.WorkingSetSize/(1<<20));
}
#endif

extern DWORD color[];

// random number generator - like drand48() under Linux, but faster
// works only with compilers which have long long int!
long long unsigned int xdr48=0x000100010001LL, mndr48=0x0005deece66dLL, doddr48=0xbLL ;
double _drand48(void)  
{
  xdr48=mndr48*xdr48+doddr48 ; xdr48&=0xffffffffffffLL ;
  return (xdr48/281474976710656.0) ;
}

void _srand48(unsigned int x0)
{
	xdr48=(x0<<16)+0x330e ; 
}

double gauss48(){

  static int i=0;
  double phi,r;
  static double x,y;

  if(i==0){
    phi = 2.0*M_PI*_drand48();
    r = sqrt(-2.0*log(_drand48()));
    x = cos(phi)*r;
    y = sin(phi)*r;
    i=1;
    return x;
  }
  else{
    i=0;
    return y;
  }
} 

double sign(double x) { if (x<0) return -1 ; else return 1 ; }

double carpol(double x,double y)
{
  if (x==0) if (y>0) return M_PI*0.5 ; else return 1.5*M_PI ;
  double f=atan(y/x);
  if (x<0) f+=M_PI ;
  if (x>0 && y<0) f+=2*M_PI ;
  return f ;
}

double ** duplicate(double **a, int x, int y)
{
  double **b=new double*[x] ; 
  for (int i=0;i<x;i++) { 
    b[i]=new double[y] ;
    for (int j=0;j<y;j++) b[i][j]=a[i][j] ;
  }
  return b ;
}

double ** allocate(int x, int y, double fill)
{
  double **b=new double*[x] ; 
  for (int i=0;i<x;i++) { 
    b[i]=new double[y] ;
    for (int j=0;j<y;j++) b[i][j]=fill ;
  }
  return b ;
}

int ** allocate(int x, int y, int fill)
{
  int **b=new int*[x] ; 
  for (int i=0;i<x;i++) { 
    b[i]=new int[y] ;
    for (int j=0;j<y;j++) b[i][j]=fill ;
  }
  return b ;
}

double ** allocate()
{
  double **b=new double*[tabx] ; 
  for (int i=0;i<tabx;i++) { 
    b[i]=new double[taby] ;
  }
  return b ;
}


void remove(double **a,int x, int y)
{
  for (int i=0;i<x;i++) delete [] a[i] ;
  delete [] a ;  
}

void deallocate_tabs() 
{
  if (cc==NULL) err("cc=null") ;
  remove(cc,tabx,taby) ; remove(dcc,tabx,taby) ; remove(sinkcc,tabx,taby) ; remove(slime,tabx,taby) ;
  delete [] ccone ; delete [] dccone ;  
}

void reallocate_tabs()
{
  int i,j,k;

  double **cc0, **slime0, **sinkcc0, *ccone0 ;
  cc0=cc ; slime0=slime ; sinkcc0=sinkcc ; ccone0=ccone ;
  remove(dcc,tabx,taby) ; delete [] dccone ;

  tabx*=2 ; taby*=2 ; heightone*=2 ; 
  cc=allocate() ; dcc=allocate() ; sinkcc=allocate() ; slime=allocate() ; 
  ccone=new double[heightone] ; dccone=new double[heightone] ;
  for (i=0;i<heightone;i++) ccone[i]=c0 ;  // no copying from old table - works fine if heightone=taby at all times
// BW heightone and ccone not updated!
  
  for (i=0;i<tabx;i++) {
    for (j=0;j<taby;j++) { 
      if (i>=tabx/4 && i<3*tabx/4 && j>=taby/4 && j<3*taby/4) {
        slime[i][j]=slime0[i-tabx/4][j-taby/4] ;        
        cc[i][j]=cc0[i-tabx/4][j-taby/4] ; sinkcc[i][j]=sinkcc0[i-tabx/4][j-taby/4] ;
      } else {
        slime[i][j]=0 ; sinkcc[i][j]=0 ;
        float r2=SQR(i+0.5-tabx/2.)+SQR(j+0.5-tabx/2.) ;
        if (int(sqrt(r2)-tabx/4)>=heightone/2) err("<0") ;
        if (r2<heightone*heightone/4) cc[i][j]=ccone0[int(sqrt(r2)-tabx/4)] ;
        else cc[i][j]=c0 ; 
      }      
    }
  }
  
  remove(cc0,tabx/2,taby/2) ; remove(sinkcc0,tabx/2,taby/2) ; remove(slime0,tabx/2,taby/2) ; delete [] ccone0; 
}

void allocate_tabs()
{
  int i,j,k;
  
  cc=allocate() ; dcc=allocate() ; sinkcc=allocate() ; slime=allocate() ; 
  ccone=new double[heightone] ; dccone=new double[heightone] ;
  for (i=0;i<heightone;i++) ccone[i]=c0 ;  

#ifdef RECTBOX
  bulkbact=new int *[height] ;
  for (j=0;j<height;j++) {
    bulkbact[j]=new int[width] ;
    for (i=0;i<width;i++)
      bulkbact[j][i]=0 ;     
  }
#endif

  for (i=0;i<tabx;i++) {
    for (j=0;j<taby;j++) { 
      slime[i][j]=0 ; sinkcc[i][j]=0 ;
#ifdef RECTBOX
      cc[i][j]=c0 ; 
#elif defined EXPANDING_COLONY
      cc[i][j]=c0 ; //*(1-0.7*exp(-0.001*(SQR(i-tabx*0.5)+SQR(j-tabx*0.5)))) ; 
#endif
    }
  }
}

inline int xtobox(double x) { return int(nbox/2+x/boxwidth) ; }
inline int ytobox(double y) { return int(nbox/2+y/boxwidth) ; }

void print_boxes()
{
	for (int i=0;i<nbox;i++) {
		for (int j=0;j<nbox;j++) cout<<boxes[i][j].size()<<" " ;
		cout<<endl ;
	}
}

void   populate_boxes()
{
  int i,j;
  for (i=0;i<nbox;i++) for (j=0;j<nbox;j++) boxes[i][j].clear() ;
  for (i=0;i<bact.size();i++) {
    bact[i]->bi=boxes[xtobox(bact[i]->x())][ytobox(bact[i]->y())].size() ;
    boxes[xtobox(bact[i]->x())][ytobox(bact[i]->y())].push_back(i) ;
    bact[i]->bx=xtobox(bact[i]->x()) ;
    bact[i]->by=ytobox(bact[i]->y()) ;
    if (bact[i]->bx<0 || bact[i]->by<0 || bact[i]->bx>=nbox || bact[i]->by>=nbox) err("bxby outside nobx") ;
  }
}

void update_box(Bacterium *b, int i)
{
  int bx=xtobox(b->x()) ;
  int by=ytobox(b->y()) ;
  if (bx<0 || by<0 || bx>=nbox || by>=nbox) err("bxby outside nobx") ;
  if (b->bi==-1) { // freshly made bacterium
    //err("bx or by=-1",i) ;
    b->bx=bx ; b->by=by ;
    b->bi=boxes[bx][by].size() ;
    boxes[bx][by].push_back(i) ;
//    return ;
  } else if (b->bx!=bx || b->by!=by) {
    if (boxes[b->bx][b->by].size()==0) err("zero",b->by) ;
    int j=boxes[b->bx][b->by].back() ;  // index of bacterium at the end of the box
    if (j<0 || j>=bact.size()) err("j=",j) ;
    //err("x",j) ;
    boxes[b->bx][b->by].pop_back() ;
    if (j!=i) { boxes[b->bx][b->by][b->bi]=j ; bact[j]->bi=b->bi ; }
    
    b->bi=boxes[bx][by].size() ; 
    b->bx=bx ; b->by=by ;
    boxes[bx][by].push_back(i) ;
  } 
}


void just_mutated(int i) 
{
  float yp=bact[i]->y(),ym=yp ;
  for (int j=0;j<bact.size();j++) {
    if (i!=j && fabs(bact[i]->x()-bact[j]->x())<2*d0 && bact[j]->y()>yp) yp=bact[j]->y() ;
    if (i!=j && fabs(bact[i]->x()-bact[j]->x())<2*d0 && bact[j]->y()<ym) ym=bact[j]->y() ;
  }
  char txt[256] ;
  sprintf(txt,"%s/runs_%d.dat",NUM,RAND) ; runs.open(txt,ios::app) ; 
  runs << "mutated "<<i<<" at_position "<<bact[i]->x()<<" "<<bact[i]->y()<<" y+ "<<yp<<" y- "<<ym<<" " ;
  runs.close() ;  
}

int Bacterium::tot_number_of_bacteria = 0 ; // initialize tot. bacterium counter

void reset(bool reset_rand)
{
	//printf("reset\n") ;
  dt=dtmax ;
  int i,j,k;
  char txt[256] ;
  for (i=0;i<bact.size();i++) delete bact[i] ;
  bact.clear() ;
  if (cc!=NULL) deallocate_tabs() ;
#ifdef EXPANDING_COLONY
  width=height=100 ;
#endif  
  dx0=sqrt(diff_rate*dt/0.125*0.5) ; if (dx0<d0) dx0=d0 ;
  tabx=int(width/dx0) ; taby=int(height/dx0) ; dx0=1.*width/tabx ; //err("x",dx0) ;
  heightone=taby ; 
    
  allocate_tabs() ;     
  
#ifdef RECTBOX  

//  float xc=-width/2+2*d0,yc=-10*d0 ;
#if defined(ECOLI)
	for (int j=0;j<0.2*width*height;j++) {
    Rod_shaped *b=new Rod_shaped(0) ; //i%2) ; 
    bool too_close=false ;
    int iter=100 ;
		do {
			b->r.x=(_drand48()-0.5)*width ; b->r.y=(_drand48()-0.5)*height ; 
#ifdef THREE_DIM
			b->r.z=_drand48()*(depth-d0) ;
#endif
    	b->phi=M_PI*_drand48() ; //M_PI/2 ; 
    	b->recalc_cspt() ;
    	b->f.zero() ;
    	b->find_acceleration(j) ;
    	too_close=false ;
    	for (i=0;i<bact.size();i++) {
				double r2=SQR(bact[i]->x()-b->r.x)+SQR(bact[i]->y()-b->r.y) ;
#ifdef TWO_DIM
				if (r2<SQR(1.6)) { too_close=true ; break ; }
#else
				if (r2<SQR(1.6) && fabs(bact[i]->z()-b->r.z)<d0) { too_close=true ; break ; }
#endif
			}
			iter-- ; if (iter<=0) break ;
    } while (!b->f.iszero() || too_close) ;
		if (iter>0) bact.push_back(b) ;  // the only place where "Rod_shaped" appears in the main code
#else
    #error this type of cell not implemented!
#endif
  } 
#elif defined EXPANDING_COLONY

  if (INITIAL_CONDITION==1) {
// single colony
    bact.push_back(new Rod_shaped(0)) ;
  } else if (INITIAL_CONDITION==2) {
//  two colliding colonies
    float l=30 ;
    Rod_shaped *bb=new Rod_shaped(0) ;
    bb->r.x=-l/2 ; bb->type=0 ;
    bact.push_back(bb) ;  
    bb=new Rod_shaped(0) ;
    bb->r.x=l/2 ; bb->type=1 ;
    bact.push_back(bb) ;
  } else if (INITIAL_CONDITION==3) {
// ring of cells
    Rod_shaped *bb ;
    int iw=0;
    for (double a=0;a<2*M_PI;a+=0.1) {
      bb=new Rod_shaped(0);
      bb->r.x=cos(a)*15 ; bb->r.y=sin(a)*15 ; bb->phi=a ; bb->recalc_cspt() ; bb->type=iw%2 ; 
      bact.push_back(bb); 
      iw++;
    }
  } else if (INITIAL_CONDITION==4) {
    // cells in a disk, separated by some distance
    Rod_shaped *bb ; 
    const double rr=70 ;
    for (double y=-rr;y<=rr;y+=3) {
    	for (double x=-rr;x<=rr;x+=3) {
    		if (x*x+y*y<rr*rr && _drand48()<0.1) {
      		bb=new Rod_shaped(0) ;
      		bb->phi=_drand48()*2*M_PI ; bb->recalc_cspt() ; 
      
					bb->r.x=(_drand48()-0.5)*0.6+x ; 
      		bb->r.y=(_drand48()-0.5)*0.6+y ;
      		bact.push_back(bb); 
      	}
      }
    }
  } else { err("INITIAL_CONDITION not defined") ; }


#else
  #error nor RECTBOX nor EXPANDING_COLONY defined
#endif

  populate_boxes() ;
 	bact[0]->init() ; // initialize some tables inside "Rod_shaped" and find closest n.n.

  if (sectors) { sprintf(txt,"%s/sectors_%d_%d.dat",NUM,RAND,sample) ; secs.open(txt) ; secs.close() ; }
  sprintf(txt,"%s/par(t)_%d_%d.dat",NUM,RAND,sample) ; params.open(txt) ; params.close() ;
  
  sprintf(txt,"%s/runs_%d.dat",NUM,RAND) ; runs.open(txt,ios::app) ; 
  if (sample>0) runs<<endl ;
  runs<<sample<<" " ; 
  runs.close() ;

  long long unsigned int temp ;   
  if (!reset_rand) temp=xdr48 ; 
  if (start_file!=NULL) { printf("loading file %s",start_file) ; fflush(stdout) ; load_all(start_file) ; printf(" done\n") ; fflush(stdout) ; }
  if (!reset_rand) xdr48=temp ; // make sure that after resetting and reloading init. conf. we do not restart RNG
//  load_all("start_circles_100x100.dat") ;
//  for (j=0;j<bact.size();j++) bact[j]->type=drand48()*2 ;

  ttt_last_save=-(1<<30) ; ttt_last_dt_update=0 ; ttt=0 ; tt=0 ; total_y_shift=0 ; prev_front_pos=0 ;
  mut_occurred=false ;
  reassign=false ;
  no_of_mutations=0 ;
  mean_front_roughness=0 ;
#ifdef JUNCTIONS
  junctions.clear() ;
#endif
#ifdef MUTATION_MODE
  if (MUTATION_MODE==2) { // only a single cell mutates upon division. all types initially set to 0
    for (i=0;i<bact.size();i++) {
      bact[i]->type=0 ; bact[i]->fitness=1 ;
    }  
    mut_occurred=false ;
  }
  if (MUTATION_MODE==3) { // a single cell is selected at random from the first line of cells and it mutates
    for (i=0;i<bact.size();i++) { bact[i]->type=0 ; bact[i]->fitness=1 ; }
    for (int t=0;t<100;t++) {
      float x0=_drand48()*width-width/2, y0=-1e6 ;
      j=-1 ;
      for (i=0;i<bact.size();i++) {
        if (bact[i]->x()<x0+d0 && bact[i]->x()>x0-d0 && bact[i]->y()>y0) { y0=bact[i]->y() ; j=i ; } 
      }
      if (j!=-1) break ;
      printf("could not select a bacterium for mutation, trying again...\n") ;
    }
    bact[j]->type=1 ; bact[j]->fitness=1+fit_adv ; just_mutated(j) ;
    mut_occurred=true ;
  }
  if (MUTATION_MODE==4) { //all cells are assigned different types from 0 to no_of_cells-1
    for (i=0;i<bact.size();i++) bact[i]->type=i ;
    mut_occurred=true ;
  }
  if (MUTATION_MODE==5) { //cells are assigned types 0 and 1 at random
    for (i=0;i<bact.size();i++) { 
			bact[i]->type=_drand48()<0.9?0:1 ; 
			//if (bact[i]->x()-0.33*bact[i]->y()>30) bact[i]->type=1 ; // BW
			if (bact[i]->type==1) bact[i]->fitness=1+fit_adv ;
		}
    mut_occurred=true ;
  }
  if (MUTATION_MODE==6) { // a random cell is selected for mutation
    for (i=0;i<bact.size();i++) { bact[i]->type=0 ; bact[i]->fitness=1 ; }
    j=int(_drand48()*bact.size()) ;
    bact[j]->type=1 ; bact[j]->fitness=1+fit_adv ; just_mutated(j) ;
    mut_occurred=true ;  
  }
#endif

  noncriterr=0 ;
}

// initialize everything: allocate memory, set initial concentrations of food
// throw some cells at random
void init()
{
  if (int(width/boxwidth)*boxwidth!=width) err("width must be the multiple of 10*d0") ;
  dt=dtmax ;
  int i,j,k;

#ifdef RECTBOX
  nbox=2+MAX(int(width/boxwidth),int(height/boxwidth)) ;
  nbox=int(nbox/2)*2 ; // to make it even
#else
  nbox=128 ; // BW: fixed size - add expanding it later!
#endif
  boxes.resize(nbox) ;
  for (i=0;i<nbox;i++) boxes[i].resize(nbox) ;

  char txt[256] ;
  sprintf(txt,"mkdir %s",NUM) ;
  system(txt) ;

  sprintf(txt,"%s/runs_%d.dat",NUM,RAND) ; runs.open(txt) ; runs.close() ;

  reset(false) ;
  
  tinit=clock() ; noncriterr=0 ;
}

#define FW(x) fwrite(&(x),sizeof((x)),1,all) ;
#define FR(x) fread(&(x),sizeof((x)),1,all) ;

void save_all(char *name)
{
  FILE *all;
  int i;
  all=fopen(name,"wb") ;
  i=bact.size() ;
  FW(i) ; FW(tabx) ; FW(taby) ; FW(heightone) ; FW(width) ; FW(height) ;
  FW(tt) ; FW(total_y_shift) ; FW(ycm) ; FW(yfront) ; FW(xdr48) ;

  for (i=0;i<bact.size();i++) bact[i]->save_data(all) ;
  for (i=0;i<tabx;i++) fwrite(cc[i],sizeof(cc[i][0]),taby,all) ;
  for (i=0;i<tabx;i++) fwrite(sinkcc[i],sizeof(sinkcc[i][0]),taby,all) ;
#ifdef RECTBOX
  for (i=0;i<height;i++) fwrite(bulkbact[i],sizeof(bulkbact[i][0]),width,all) ;
#endif
  for (i=0;i<tabx;i++) fwrite(slime[i],sizeof(slime[i][0]),taby,all) ;
  fwrite(ccone,sizeof(double),heightone,all) ;
  fclose(all) ;
}

void load_all(char *name)
{
  int i,nn;
  int tabxn, tabyn, heightonen ;
  long long int newxdr48 ;
  FILE *all;
  all=fopen(name,"rb") ;
  FR(nn) ; FR(tabxn) ; FR(tabyn) ;FR(heightonen) ; FR(width) ; FR(height) ;
  FR(tt) ; FR(total_y_shift) ; FR(ycm) ; FR(yfront) ; FR(newxdr48) ;
  if (tabxn!=tabx || tabyn!=taby || heightonen!=heightone) { deallocate_tabs() ; tabx=tabxn ; taby=tabyn ; heightone=heightonen ; allocate_tabs() ; }
  for (i=0;i<bact.size();i++) delete bact[i] ;
  bact.clear() ;
  for (i=0;i<nn;i++) { 
    bact.push_back(new Rod_shaped) ; bact[i]->load_data(all) ;
  }
#ifdef JUNCTIONS
  junctions.clear() ; // BW this should be replaced with loading junctions from saved file but first I would have to replace pointers with indices for b1,b2
#endif
  for (i=0;i<tabx;i++) fread(cc[i],sizeof(cc[i][0]),taby,all) ;
  for (i=0;i<tabx;i++) fread(sinkcc[i],sizeof(sinkcc[i][0]),taby,all) ;
#ifdef RECTBOX
  for (i=0;i<height;i++) fread(bulkbact[i],sizeof(bulkbact[i][0]),width,all) ;
#endif
  for (i=0;i<tabx;i++) fread(slime[i],sizeof(slime[i][0]),taby,all) ;
  fread(ccone,sizeof(double),heightone,all) ;
      
  fclose(all) ; 
  populate_boxes() ;
  xdr48=newxdr48 ;
  _maxdisp=0 ;
}

void save_data(char *name, int cc_all=0)
{
  FILE *pos=fopen(name,"w") ;
  if (pos==NULL) err(name) ;
  fprintf(pos,"%d %d %d\t%d %d\n",bact.size(),tabx,taby,width,height) ;

  for (int i=0;i<bact.size();i++) {
    char txt[256] ;
    if ((save_cells&64)) bact[i]->to_string_avvel(txt) ; else bact[i]->to_string(txt) ;
    fprintf(pos,txt) ;
  }
  if (cc_all==0) {
    for (int i=0;i<tabx;i++) {
      for (int j=0;j<taby;j++) 
        fprintf(pos,"%lf ",cc[i][j]) ;
      fprintf(pos,"\n") ;
    }  
  }
  fclose(pos) ;
}


double globalOrder() { //output the global order parameter

     double Q11=0,Q12=0;

     for (int i=0;i<bact.size();i++) {
         Q11+=SQR(bact[i]->cos_phi())-0.5;
         Q12+=bact[i]->cos_phi() * bact[i]->sin_phi();
     }

     Q11 /= bact.size(); Q12 /= bact.size();

     return 2*sqrt(SQR(Q11) + SQR(Q12));
}

int count_types(int t)
{
  int n=0 ;
  for (int i=0;i<bact.size();i++) if (bact[i]->type==t) n++ ;
  return n ;
}

bool is_only_one_type()
{
  for (int i=1;i<bact.size();i++) if (bact[i]->type!=bact[0]->type) return false ;
  return true ; 
}

void save_how_many_bact_each_type(char *name)
{
  int maxt=0 ;
  for (int i=0;i<bact.size();i++) if (bact[i]->type+1>maxt) maxt=bact[i]->type+1 ;
  int *types=new int[maxt], *act_types=new int[maxt] ;
  for (int i=0;i<maxt;i++) types[i]=act_types[i]=0 ;
  for (int i=0;i<bact.size();i++) {
    types[bact[i]->type]++ ;
    if (bact[i]->active==1) act_types[bact[i]->type]++ ;
  }

  FILE *f ;
  if (ttt==0) f=fopen(name,"w") ; else f=fopen(name,"a") ;
  if (f==NULL) err("cannot open file for save_how_many_bact_each_type") ;
  for (int i=0;i<maxt;i++) if (types[i]>0) fprintf(f,"%lf %d %d %d\n",tt,i,types[i],act_types[i]) ; // time_step  type  #bacteria
  fclose(f) ;
  delete [] types ; delete [] act_types ;
}

void save_positions(char *name)
{

  FILE *f ;
  if (ttt==0) f=fopen(name,"w") ; else f=fopen(name,"a") ;
  if (f==NULL) err("cannot open file for save_positions") ;
  char txt[256] ;
  for (int i=0;i<bact.size();i++) {
    bact[i]->to_string_avvel(txt) ; fprintf(f,"%lf %lf %d %d\t",tt,prev_front_pos,bact[i]->number,bact[i]->prev_number) ; fprintf(f,txt) ;
  }
  fclose(f) ;
}


#ifdef TUBE
float roughness()
{
  int w=int(width),j ;
  float *h=new float[w] ;
  for (int i=0;i<w;i++) h[i]=2e6 ; // something big
  for (int i=0;i<bact.size();i++) {
    vec r ;
    bact[i]->fill_lowest_y(h,w) ;
    //bact[i]->lowest_y(&r) ; j=int(r.x+w/2+w)%w ; if (r.y<h[j]) h[j]=r.y ;
  } 
  float sum=0,sum2=0 ;
  int nn=0 ;
  for (int i=0;i<w;i++) if (h[i]<1e6) { sum+=h[i] ; sum2+=h[i]*h[i] ; nn++ ; } 
  delete [] h; 
  if (nn==0) return 0 ; else { mean_front_roughness=sqrt((sum2-sum*sum/nn)/nn) ; return mean_front_roughness ; }
}

float * map_front_thickness()  // returns a pointer to an array of thickness
{
  min_front_thickness=1e6 ; mean_front_thickness=0 ;
  int w=int(width),j ;
  float *hh=new float[w],*hm=new float[w] ;
  for (int i=0;i<w;i++) { hm[i]=2e6 ; hh[i]=-2e6 ; } // something big
  for (int i=0;i<bact.size();i++) {
    vec r ;
    bact[i]->fill_lowest_y(hm,w) ;
    bact[i]->fill_highest_y(hh,w) ;
//    bact[i]->lowest_y(&r) ; j=int(r.x+w/2+w)%w ; if (r.y<hm[j]) hm[j]=r.y ;
//    bact[i]->highest_y(&r) ; j=int(r.x+w/2+w)%w ; if (r.y>hh[j]) hh[j]=r.y ;
  } 
  for (int i=0;i<w;i++) { hh[i]-=hm[i] ; mean_front_thickness+=hh[i] ; if (min_front_thickness>hh[i]) min_front_thickness=hh[i] ; }
  mean_front_thickness/=w ;
  delete [] hm ;   
  return hh ;
}

float thickness()
{
  int w=int(width),j ;
  float *hh=new float[w],*hm=new float[w] ;
  for (int i=0;i<w;i++) { hm[i]=2e6 ; hh[i]=-2e6 ; } // something big
  for (int i=0;i<bact.size();i++) {
    vec r ;
    if (bact[i]->active>0) {
      bact[i]->fill_lowest_y(hm,w) ;
      bact[i]->fill_highest_y(hh,w) ;
    }
//    bact[i]->lowest_y(&r) ; j=int(r.x+w/2+w)%w ; if (r.y<hm[j]) hm[j]=r.y ;
//    bact[i]->highest_y(&r) ; j=int(r.x+w/2+w)%w ; if (r.y>hh[j]) hh[j]=r.y ;
  } 
  float dh=0 ;
  int nn=0 ;
//  for (int i=0;i<w;i++) if (hh[i]!=-2e6 && hm[i]!=2e6) { dh+=hh[i]-hm[i] ; nn++ ; }
  for (int i=0;i<w;i++) if (hh[i]!=-2e6) {
    float dhmin=1e6 ;
    for (int j=-0.1*w;j<0.1*w;j++) if (hm[(w+j+i)%w]!=2e6) {
      float d=sqrt(SQR(1.*j)+SQR(hh[i]-hm[(w+j+i)%w])) ;
      if (d<dhmin) dhmin=d ;
    } 
    dh+=dhmin ; nn++ ;
  }
/*  FILE *f=fopen("test.dat","w") ;
  for (int i=0;i<w;i++) fprintf(f,"%d %f\n",i,hh[i]) ;
  for (int i=0;i<w;i++) fprintf(f,"%d %f\n",i,hm[i]) ;
  fclose(f) ;*/
  delete [] hh; delete [] hm ;
//  if (dh<0) dh=0 ; 
  return dh/nn ;
}

double filling_fraction()
{
  int occup=0 ;
  for (int i=1;i<=10;i++) for (int j=0;j<width;j++) if (bulkbact[height-i][j]>0) occup++ ;
  return 1.*occup/(10*tabx);
}

#endif

#ifdef EXPANDING_COLONY
float radius()
{
  float rr[10],rad=0 ;
  for (int i=0;i<10;i++) rr[i]=-1 ;
  for (int i=0;i<bact.size();i++) {
    float phi=carpol(bact[i]->x(),bact[i]->y()) ;
    float r=sqrt(SQR(bact[i]->x())+SQR(bact[i]->y())) ;
    while (phi<0) phi+=2*M_PI ; while (phi>2*M_PI) phi-=2*M_PI ;
    int j=int(10*phi/(2*M_PI)) ;
    if (rr[j]<r) rr[j]=r ;
  }
  int n=0 ;
  for (int i=0;i<10;i++) if (rr[i]>0) { rad+=rr[i] ; n++ ; }
  if (n>0) return rad/n ;
  return 0 ;
}
#endif

int active_cells()
{
  int nn0=0 ;
  for (int i=0;i<bact.size();i++) if (bact[i]->active==1) nn0++ ;
  return nn0 ;
}  

int growing_cells()
{
  float max_gr_rate=0 ;
  int nn=0 ;
  for (int i=0;i<bact.size();i++) if (bact[i]->elong_rate()>max_gr_rate) max_gr_rate=bact[i]->elong_rate() ;
  for (int i=0;i<bact.size();i++) if (bact[i]->elong_rate()>0.5*max_gr_rate) nn++ ;
  return nn ;
}

inline int xtotab(double x) { return int(tabx/2+x/dx0) ; }
inline int ytotab(double y) { return int(taby/2+y/dx0) ; }


void save_profile(char *name)
{
  int h=int(height/d0) ;
  float *profile=new float[h] ; // concentration of cells along y-axis
  int i,j;
  for (j=0;j<h;j++) profile[j]=0 ; 
  
  for (i=0;i<bact.size();i++) {
    j=int(h/2+bact[i]->y()/d0) ;
    profile[j]+=1./width ; 
  }
//  static int first=1 ;
  FILE *prof ;
//  if (first) {
  prof=fopen(name,"w") ; //first=0 ;
//  } else prof=fopen(name,"a") ;
  for (int j=0;j<h;j++) fprintf(prof,"%.3f ",profile[j]) ; fprintf(prof,"\n") ;
  fclose(prof) ;
  delete [] profile ;
}

void save_envelope(char *name)
{
  int w=int(width/d0) ;
  float *env=new float[w] ; // envelope of the front
  int i,j ;
  for (i=0;i<w;i++) env[i]=height ;
//  static int first=1 ;
  for (i=0;i<bact.size();i++) {
    j=int(w/2+bact[i]->x()/d0) ;
    if (bact[i]->y()<env[j]) env[j]=bact[i]->y() ; 
  }
  FILE *en ;
//  if (first) {
  en=fopen(name,"w") ; //first=0 ;
//  } else en=fopen(name,"a") ;
  for (i=0;i<w;i++) fprintf(en,"%.3f ",-env[i]) ; fprintf(en,"\n") ;
  fclose(en) ;
  delete [] env ;
}

void end()
{
// if any files open, close them
}


//inline int xtotab(double x) { int i=int(3*tabx/2+x/dx0)%tabx ; /*if (i<0) err("i<0",i) ; if (i>=tabx) err("i>=tabx",i) ;*/ return i ; }
//inline int ytotab(double y) { int j=int(3*taby/2+y/dx0)%taby ; /*if (j<0) err("j<0",j) ; if (j>=taby) err("j>=taby",j) ;*/ return j ; }


#ifdef NUTRIENT
// solves the diffusion equation numerically
void diffusion()
{
  int i,j,i0,i1,j0,j1;
  i0=0 ; i1=tabx ; j0=0; j1=taby ;
  for (i=1;i<tabx-1;i++)
    for (j=1;j<taby-1;j++) 
      dcc[i][j]=diff_rate*(cc[i+1][j]+cc[i-1][j]+cc[i][j+1]+cc[i][j-1] -4*cc[i][j])/(dx0*dx0);
  // bottom row - no flux condition
  j=taby-1 ; for (i=1;i<tabx-1;i++) dcc[i][j]=diff_rate*(cc[i+1][j]+cc[i-1][j]+cc[i][j-1] -3*cc[i][j])/(dx0*dx0);
  // left and right - no flux condition
  i=0 ; for (j=0;j<taby-1;j++) dcc[i][j]=diff_rate*(cc[i+1][j]+cc[i][j+1]+cc[i][j-1] -3*cc[i][j])/(dx0*dx0);
  i=tabx-1 ; for (j=0;j<taby-1;j++) dcc[i][j]=diff_rate*(cc[i-1][j]+cc[i][j+1]+cc[i][j-1] -3*cc[i][j])/(dx0*dx0);
	// first upper row
  for (i=0;i<tabx;i++) dcc[i][0]=diff_rate*(c0+cc[i][1] -2*cc[i][0])/(dx0*dx0);
	// corners
	i=0 ; j=taby-1 ; dcc[i][j]=diff_rate*(cc[i+1][j]+cc[i][j-1] -2*cc[i][j])/(dx0*dx0);
	i=tabx-1 ; j=taby-1 ; dcc[i][j]=diff_rate*(cc[i-1][j]+cc[i][j-1] -2*cc[i][j])/(dx0*dx0);

  for (i=i0;i<i1;i++)
    for (j=j0;j<j1;j++) {
      cc[i][j]+=dt*dcc[i][j] ; 
      if (sinkcc[i][j]>0) { 
        cc[i][j]-=dt*sinkcc[i][j]*uptakefun(cc[i][j]) ;
      }
      if (cc[i][j]<-c0 || cc[i][j]>20*c0) noncriterr=1 ;
      if (cc[i][j]<0) cc[i][j]=0 ;
    }
//  for (i=0;i<tabn;i++) cc[i][tabn/2]=0 ;
}
#endif

void rescale() ;
//==============================================================================

// main loop - does a single step of time evolution
void run()
{
  if (_stop || bact.size()==0 || bact.size()>max_size) return ;


  int i,j,k;

  how_many_steps_save=int(1/dt) ; 

#ifdef TUBE
  if (ttt%10==0) {
    ycm=0 ; // center of mass
    for (i=0;i<bact.size();i++) { ycm+=bact[i]->y()/bact.size() ; } //if (bact[i]->y()<y0) y0=bact[i]->y() ; }
    yfront=1e6 ; k=0; // yfront is the minimal (highest) value of y at the front 
    for (i=0;i<bact.size();i++) if (bact[i]->y()<yfront) yfront=bact[i]->y() ;
  }
#endif
    
  _maxdisp=0 ; nointeracting=0 ;
#ifdef CLOSEST  
  noclosest=0 ;
#endif
  
  for (i=0;i<bact.size();i++) {
    bact[i]->recalc_mass_and_I() ;
    bact[i]->clear_forces() ; 
  }
  for (i=0;i<bact.size();i++) bact[i]->find_acceleration(i) ;  
#ifdef JUNCTIONS
  simulate_junctions() ; 
#endif

#ifdef OVERDAMPED
  for (i=0;i<bact.size();i++) if (!noncriterr && bact[i]->active>=0) bact[i]->one_Euler_step(i) ; 
#else
  for (i=0;i<bact.size();i++) if (!noncriterr && bact[i]->active>=0) bact[i]->one_Newton_step(i) ; 
#endif

//  for (i=0;i<bact.size();i++) if (isinf(bact[i]->x()) || isinf(bact[i]->y())) err("isinf") ;
  
  if (!noncriterr) { // if everything is ok then....
    int nn0=bact.size() ;
    if (gr_ok) { for (i=0;i<nn0;i++) bact[i]->one_growth_step(i) ; }
    
    for (i=0;i<bact.size();i++) update_box(bact[i],i) ;

#ifdef CLOSEST    
    avclosest=0 ;
    for (i=0;i<bact.size();i++) {
      if (bact[i]->nn_needs_update && bact[i]->active>=0) bact[i]->find_closest(i) ;
      avclosest+=bact[i]->closest.size() ;
    }
    avclosest/=bact.size() ;
#endif
#ifdef NUTRIENT
    diffusion() ; // lastly, let's do diffusion
#endif
  } else if ((save_cells&1) || backup) { // but if not, load the last configuration and decrease dt
    char txt[256] ; 
    if (backup) sprintf(txt,"%s/temp_%d.dat",NUM,RAND) ; 
    if ((save_cells&1)) sprintf(txt,"%s/temp_%d_%d.dat",NUM,RAND,sample) ; 
    printf("dt decreased at t=%lf to %lf, reloading configuration...",tt,dt/2) ; fflush(stdout) ;
    load_all(txt) ; printf(" done\n") ; fflush(stdout) ;
    dt/=2 ; noncriterr=0 ;
    if (dt<1e-6) err("dt too small!") ;
    return ;    
  } else {//err("dt too large. Cannot reload, config. non-existent") ;
    // reinitialize if nothing else can be done, with smaller dt
    printf("cannot decrease dt because temp.dat does not exists. Restarting...\n") ; fflush(stdout) ;
    reset(true) ;
    dt/=2 ; noncriterr=0 ;
    return ;
  }
  
#ifdef TWO_DIM
  if (fmag_transition>0) {
    for (i=0;i<bact.size();i++) {
      if (bact[i]->fmagav>fmag_transition) { bact[i]->active=-1 ; }  // remove the bacterium if too compressed
    }
  }
#endif

  if (ttt%100==0) {
    xmin=1e6 ; xmax=-1e6 ;ymin=1e6 ;ymax=-1e6 ;  
    zmin=1e+6 ; zmax=-1e-6 ;
    for (i=0;i<bact.size();i++) {
      Bacterium *b=bact[i];
      if (b->x()>xmax) xmax= b->x() ;
      if (b->x()<xmin) xmin= b->x() ;
      if (b->y()>ymax) ymax= b->y() ;
      if (b->y()<ymin) ymin= b->y() ;
      if (b->z()>zmax) zmax=b->z() ;
      if (b->z()<zmin) zmin=b->z() ;

    }

#ifdef CLOSEST
    int some_change ;
    some_change=0 ;
  // change state of cells which do not grow
    for (i=0;i<bact.size();i++) {
//    vely[ytotab(bact[i]->y())]+=bact[i]->vy() ; velyn[ytotab(bact[i]->y())]++ ; //tab[(tabx+xtotab(x[i]))%tabx][(taby+ytotab(y[i]))%taby][0] ;
      if (bact[i]->change_state() || bact[i]->active==-1) some_change=1 ;
    }
    
    for (i=0;i<bact.size();i++) {
      if (ytotab(bact[i]->y())<=4*d0) { some_change=1 ; bact[i]->active=-1 ; }  // in addition, remove the bacterium if outside the box
//      if (ytotab(bact[i]->y())>=taby*0.65) { some_change=1 ; bact[i]->active=-1 ; }  //  BW
      if (bact[i]->active==-1) { // delete bacterium if far from the front

        if (ytotab(bact[i]->y())<0 || ytotab(bact[i]->y())>=taby) err("sinkcc too small, increase height") ;
#ifdef NUTRIENT
          sinkcc[xtotab(bact[i]->x())][ytotab(bact[i]->y())]=bact[i]->uptake/(dx0*dx0) ;
#endif
        bact[i]->fill_table(bulkbact,width,height) ;
        // first save its position etc. into the file.......
        if (sectors) { 
          char txt[256], name[256] ;
					#ifdef TUBE
          	bact[i]->add_y(total_y_shift) ;
          	bact[i]->to_string(txt) ;
          	sprintf(name,"%s/sectors_%d_%d.dat",NUM,RAND,sample) ;
          	secs.open(name,ios::app) ;
          	secs<<txt ;
          	secs.close() ;
          #elif defined RECTBOX
          	bact[i]->to_string(txt) ;
          	sprintf(name,"%s/sectors_%d_%d.dat",NUM,RAND,sample) ;
          	secs.open(name,ios::app) ;
          	secs<<tt<<" "<<txt ;
          	secs.close() ;          
          #endif
        }

        // .......and then delete it:

//        if (bact[i]->bi==-1 ||  bact[i]->bx==-1 || bact[i]->by==-1) err("2. bixy==-1") ;
        j=boxes[bact[i]->bx][bact[i]->by].back() ; 
        if (j<0 || j>bact.size()) err("er. j=",j) ;
        boxes[bact[i]->bx][bact[i]->by].pop_back() ;
        if (j!=i) { boxes[bact[i]->bx][bact[i]->by][bact[i]->bi]=j ; bact[j]->bi=bact[i]->bi ; }

  
        delete bact[i] ;
        if (i<bact.size()-1) {
          bact[i]=bact.back() ; 
//          if (bact[i]->bi==-1) err("3. bi==-1, i=",i) ;
          boxes[bact[i]->bx][bact[i]->by][bact[i]->bi]=i ;           
        }
        bact.pop_back() ;
        i-- ; some_change=1 ;
      }
      if (i>=bact.size()-1) break ;
    }
// BW this is necessary only if bacteria are removed, but can slow down the program by 50% (?)
    if (some_change) for (i=0;i<bact.size();i++) if (bact[i]->active>=0) bact[i]->find_closest(i) ; // update table of n.n. after some bacteria have been deleted 
//#endif  
//#ifdef TUBE
//    while (yfront<-dx0+shift_cm_y) { shift_everything_down() ; yfront+=dx0 ; } 
#endif

  }

  char txt[256];
#ifdef TUBE
  if ((save_cells&32) && fabs(prev_front_pos-(yfront+total_y_shift))>d0) { 
    prev_front_pos=(yfront+total_y_shift) ;
    sprintf(txt,"%s/positions_%d_%d.dat",NUM,RAND,sample) ; save_positions(txt) ; 
  }
#else
/*  if ((save_cells&32)) {
    float rad=radius() ;
    if (fabs(prev_front_pos-rad)>d0) { 
      prev_front_pos=rad ;
      sprintf(txt,"%s/positions_%d_%d.dat",NUM,RAND,sample) ; save_positions(txt) ; 
    }
  }*/
#endif

  if (ttt>ttt_last_save+how_many_steps_save || bact.size()>=max_size) {
    ttt_last_save=ttt ; 

    if ((save_cells&1)) { sprintf(txt,"%s/temp_%d_%d.dat",NUM,RAND,sample) ; save_all(txt) ; }
    else if (backup) { sprintf(txt,"%s/temp_%d.dat",NUM,RAND) ; save_all(txt) ; }
//    load_all(txt) ;
//    if ((save_cells&2)) { sprintf(txt,"%s/data_%d_%d_%d.dat",NUM,RAND,sample,ttt) ; save_data(txt) ; }
    if ((save_cells&2)) { sprintf(txt,"%s/data_%d_%d.dat",NUM,RAND,sample) ; save_data(txt) ; }

//    if (ttt>=40000) exit(0) ;
#ifdef TUBE
    if ((save_cells&4)) { sprintf(txt,"%s/profile_%d_%d.dat",NUM,RAND,sample) ; save_profile(txt) ; }
//    sprintf(txt,"%s/vy_%d.dat",NUM,sample) ; save_velocity(txt) ;
    if ((save_cells&8)) { sprintf(txt,"%s/envelope_%d_%d.dat",NUM,RAND,sample) ; save_envelope(txt) ; }
#endif

    if ((save_cells&16)) { sprintf(txt,"%s/types_%d_%d.dat",NUM,RAND,sample) ; save_how_many_bact_each_type(txt) ; }

#ifdef RECTBOX
  if ((save_cells&32)) { sprintf(txt,"%s/positions_%d_%d.dat",NUM,RAND,sample) ; save_positions(txt) ;  }
#endif

/*    char name[256] ;
    static int no=0 ;
    sprintf(name,"movie/data%d.dat",no++) ;
    save_data(name,1) ;
*/


    sprintf(txt,"%s/par(t)_%d_%d.dat",NUM,RAND,sample) ; params.open(txt,ios::app) ;
    // 1.time   2.iterations  3.size  4.y_shift/radius   5.order_param   6.zmax
    // 7.thickness  8.active_cells 9.cells>0.5max_gr_rate 10.roughness 11.no.of_type_0  12. filling_fraction
#ifdef RECTBOX
    params<<tt<<"\t"<<ttt<<"\t"<<bact.size()<<"\t"<<0<<"\t"<<globalOrder()<<"\t"<<zmax ;
    params<<"\t"<<0/*thickness()*/<<"\t"<<active_cells()<<"\t"<<growing_cells()<<"\t"<<0/*roughness()*/<<"\t"<<count_types(0)<<"\t"<<0/*filling_fraction()*/ ;
#else
    params<<tt<<"\t"<<ttt<<"\t"<<bact.size()<<"\t"<<radius()<<"\t"<<globalOrder()<<"\t"<<zmax ;
    params<<"\t"<<"n.a."<<"\t"<<active_cells()<<"\t"<<growing_cells()<<"\t"<<"n.a."<<"\t"<<count_types(0)<<"\t"<<"n.a." ;
#endif
    // heightone  time[s]  memory_taken[MB]
    params<<"\t"<<heightone<<"\t"<<(clock()-tinit)/CLOCKS_PER_SEC<<"\t"<<memory_taken()<<endl ;
    params.close() ;
    if (params.fail()) params.clear() ;
  }

  ttt++ ;
  tt+=dt ;

#ifdef REASSIGN_TYPES_N
  //if (bact.size()<REASSIGN_TYPES_N) reassign=false ;
  if (bact.size()>=REASSIGN_TYPES_N && reassign==false) {
    reassign=true ; 
    for (i=0;i<bact.size();i++) {
			// BW bact[i]->type=i ;
			bact[i]->type=_drand48()<0.8?0:1 ;
			if (bact[i]->type==1) bact[j]->fitness=1+fit_adv ;

			
			// this assigns colours based on y position at the time of assignment
			//float x=(bact[i]->y()+height/2)/height ;
			//color[i]=int(x*255)+(int(4095*x*x*(1-x)*(1-x))<<8)+(int(255*(1-x))<<16) ;
		}
    mut_occurred=true ;
    rescale() ;
  }
#endif

#ifdef REASSIGN_TYPES_T
  if (tt<REASSIGN_TYPES_T) reassign=false ;
  if (tt>=REASSIGN_TYPES_T && reassign==false) {
    reassign=true ; 
    for (i=0;i<bact.size();i++) {
			bact[i]->type=i ;
			// this assigns colours based on y position at the time of assignment
			//float x=(bact[i]->y()+height/2)/height ;
			//color[i]=int(x*255)+(int(4095*x*x*(1-x)*(1-x))<<8)+(int(255*(1-x))<<16) ;
		}
    mut_occurred=true ;
    rescale() ;
  }
#endif


// BW -- this inserts a mutant in the first line at a given time
/*  if (mut_occurred && tt>16) {
//    for (i=0;i<bact.size();i++) { bact[i]->type=0 ; bact[i]->fitness=1 ; }
    j=-1 ; float y0=1e6 ;
    for (i=0;i<bact.size();i++) {
      if (bact[i]->y()<y0) { y0=bact[i]->y() ; j=i ; }
    }
    bact[j]->type=1 ; bact[j]->fitness=1.05 ; //just_mutated(j) ;
    mut_occurred=false ;
  }
*/
  
#ifdef TRY_INCREASE_DT
  if (ttt>ttt_last_dt_update+how_many_steps_save*10) { 
    if (dt<dtmax) dt*=2 ; 
    ttt_last_dt_update=ttt ; 
  }
#endif  

}
