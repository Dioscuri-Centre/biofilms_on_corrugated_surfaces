#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;
#include "functions.h"
#include "vecs.h"
#include "classes.h"
#include "params.h"

//#define FRICTION
const float fr_coeff=0.5;
//#define SLIME
#ifdef NUTRIENT
  #define NUTRIENT_LIMITED
#endif
//#define WASTE_LIMITED


#ifdef SHORT_RAND_RODS
  #define RANDOMIZE_DIRS	// if defined, random directions of bacteria are chosen after each division. Makes sense only when max_length=0.5
  #define MAX_ELONG_RATE // if defined, growth_rate specifies max. elongation rate in [um/h]
#endif

// BW define this also for E coli
#define MAX_ELONG_RATE // if defined, growth_rate specifies max. elongation rate in [um/h]

const int ttt_junctions=1000 ;

#ifdef OBSTACLES
extern int **obst ; 
#endif

extern vector <Rod_shaped*> bact ;
//extern vector <Lines> lines ;
extern vector <vector <vector <int> > > boxes ;
#ifdef JUNCTIONS
extern vector <Junction> junctions ;
#endif
#ifdef AGAROSE
extern Lattice *latt ;
#endif

extern int nointeracting ;
extern float conjugation_rate ;
//inline int xtotab(double x) { return int(3*tabx/2+x/dx0)%tabx ; }
//inline int ytotab(double y) { return int(3*taby/2+y/dx0)%taby ; }

inline int xtotab(double x) { return int(tabx/2+x/dx0) ; }
inline int ytotab(double y) { return int(taby/2+y/dx0) ; }

int Rod_shaped::powmax15=1000, Rod_shaped::powmax05=1000 ;
int Rod_shaped::rods=0 ;
double *Rod_shaped::pow15=new double[Rod_shaped::powmax15+1], *Rod_shaped::pow05=new double[Rod_shaped::powmax05+1] ;

#ifdef TWO_DIM
typedef vec2 vec ;
#elif defined THREE_DIM
typedef vec3 vec ;
#else
  #error dimensionality not defined
#endif

void Rod_shaped::init()
{
  int i;
  for (i=0;i<=powmax15;i++) pow15[i]=(4./3)*E_bact*sqrt(r0*d0)*d0*pow(1.*i/powmax15,1.5) ; ; // force give by pow15 is in Pa um^2 = pN 
  for (i=0;i<=powmax05;i++) pow05[i]=pow(1.*i/powmax05,0.5) ;  
#ifdef CLOSEST
  for (int i=0;i<bact.size();i++) bact[i]->find_closest(i) ; 
#endif
}

bool Rod_shaped::inside_bact(vec p) // returns true if point p lies inside the bacterium
{
  double dd,t ;
  find_sqr_distance_to_rod(p,dd,t) ;
  if (dd<r0*r0) return true ;
  return false ;
}


/*void per_diff(vec &r1, vec &r2, vec &dr) // calculates the different between r1 and r2 assuming pbc in x
{
  dr=r1-r2 ; 
  if (dr.x<-width/2) dr.x+=width ; 
  else if (dr.x>width/2) dr.x-=width ; 
}*/

/*void per_diff(vecd &r1, vecd &r2, vecd &dr) // calculates the different between r1 and r2 assuming pbc in x
{
  dr=r1-r2 ; 
  if (dr.x<-width/2) dr.x+=width ; 
  else if (dr.x>width/2) dr.x-=width ; 
}*/



const int kx[9]={0,1,1,0,-1,-1,-1,0,1},ky[9]={0,0,1,1,1,0,-1,-1,-1} ;  
inline int xtobox(double x) { return int(nbox/2+x/boxwidth) ; }
inline int ytobox(double y) { return int(nbox/2+y/boxwidth) ; }

int xxx ;
#define SWAP(x, y) xxx=(x) ; (x)=(y) ; (y)=xxx 

#ifdef CLOSEST
  void Rod_shaped::update_old_r() 
  {
    rold1=r+n*d ; rold2=r-n*d ; //dold=d ;    
  }

  void Rod_shaped::find_closest(int thisi) 
  {
    noclosest++ ;
    vec dr ;
    rold1=r+n*d ; rold2=r-n*d ; dold=d ;
    closest.clear() ;
    
  //  for (int i=0;i<bact.size();i++) {
    int bxn, byn, *box ;  
    if (bi==-1 || bx==-1 || by==-1) { bxn=xtobox(r.x) ; byn=ytobox(r.y) ; } else { bxn=bx ; byn=by ; }
  //  if (bxn!=63 && bxn!=64) err("x",bxn) ;
    for (int k=0;k<=8;k++) {
      box=&(boxes[(bxn+kx[k]+nbox)%nbox][(byn+ky[k]+nbox)%nbox][0]) ;
      int nmax=boxes[(bxn+kx[k]+nbox)%nbox][(byn+ky[k]+nbox)%nbox].size() ;
      for (int n=0;n<nmax;n++) {
        int i=box[n] ;
      
  //    if (i!=thisi) err("x",i) ;
        dr = this->r - bact[i]->r ;
        double r2=squared(dr) ;
        if (i!=thisi && r2<(SQR(2*max_length*d0+d0))) {
          closest.push_back(i) ; //err("x",i) ;
  /*      
        if (closest.size()>1) for (int n=closest.size()-1;n>0;n--) if (closest[n]<closest[n-1]) { xxx=closest[n] ; closest[n]=closest[n-1] ; closest[n-1]=xxx ; }
  */
          int j ;
          for (j=0;j<bact[i]->closest.size();j++) if (bact[i]->closest[j]==thisi) break ;
          if (j==bact[i]->closest.size()) bact[i]->closest.push_back(thisi) ;
        }
      }  
    }
  //  for (int j=1;j<closest.size();j++) { SWAP(closest[0],closest[j]) ; }
    nn_needs_update=0 ;
  }
#endif

int is_junction(Rod_shaped *b1,Rod_shaped *b2) ;
void add_junction(Rod_shaped *b1, Rod_shaped *b2, vec p1, vec p2)  ;
void add_junction(Rod_shaped *, vec , vec )  ;


#ifdef BOTTOM_SURFACE
/*double zigzag(double x) {
	double dx=x-floor(x) ;
	if (dx<0.5) return 2*dx ;
	return 2-2*dx ;
}*/

double ysurface(double x) {
  //return 2-10*zigzag(0.35+2*x/width) ;  // 2-10
	return height/2-4*d0-bottom_ampl+bottom_ampl*cos(2*M_PI*x/bottom_period) ;
	//return height/3+0.5*height*SQR(x/width);
  //return 0*height+0.25*height*sqrt(1-SQR(2*x/width)) ;
}

inline double sqr_distance_to_line(vec &rr, vec2 &p1, vec2 &p2, double &t)
{
  double drx=rr.x-p1.x, dry=rr.y-p1.y ;
  double nx=p2.x-p1.x, ny=p2.y-p1.y ;
  t=(drx*nx+dry*ny)/(nx*nx+ny*ny) ;   
  return (SQR(drx-t*nx)+SQR(dry-t*ny)) ;
}


double find_sqr_dist_to_curve(vec r,  vec2 &p)
// calculates the approx. closest distance between point r and the curve given by ysurface
// returns the squared distance and (in p) the normal vector to the curve
{
	double t,d2 ;
	vec2 p1=vec2(r.x-0.3,ysurface(r.x-0.3)), p2=vec2(r.x+0.3,ysurface(r.x+0.3)) ;
	d2=sqr_distance_to_line(r,p1,p2,t) ;  
	//err("x",d2) ;
	p=p1+(p2-p1)*t ;
	p.x=r.x-p.x ; p.y=r.y-p.y ;
	normalize(p) ;
//	err("x",p.x) ;
	if (p.y>0) p.y=-p.y ; // this ensures the force is always away from the substrate (assuming the substrate at the bottom)
	return d2 ;
}
#endif


#ifdef THREE_DIM
inline void Rod_shaped::find_sqr_distance_to_rod(vec3 &rr, double &dist2, double &t)
{
  double drx=rr.x-r.x, dry=rr.y-r.y,drz=rr.z-r.z ;
  t=(drx*n.x+dry*n.y+drz*n.z) ;   
  if (t>d) t=d ;
  if (t<-d) t=-d ;  
  dist2=(SQR(drx-t*n.x)+SQR(dry-t*n.y)+SQR(drz-t*n.z)) ;
}

void force_from_substrate(double m, vec3 r,vec &f)
{
  f.zero() ;
	
	const double fff=1e6 ;
	if (r.x>width/2-d0/2) f.x-=fff*(r.x-(width/2-d0/2)) ;
	if (r.x<-width/2+d0/2) f.x-=fff*(r.x-(-width/2+d0/2)) ;
	if (r.y>height/2-2*d0) f.y-=fff*(r.y-(height/2-2*d0)) ;
	if (r.z>depth-3*d0/2) f.z-=fff*(r.z-(depth-3*d0/2)) ;
	if (r.z<-d0/2) f.z-=fff*(r.z-(-d0/2)) ;

#ifdef BOTTOM_SURFACE
	vec2 p;
	double d2=find_sqr_dist_to_curve(r,p) ;
	if (d2<d0*d0) { double q=fff*(sqrt(d2)-d0) ; f.x-=p.x*q ; f.y-=p.y*q ; }
	else if (r.y-d0/2>ysurface(r.x)) f.y-=fff*(r.y-d0/2-ysurface(r.x)) ;
#endif

  
}
#elif defined TWO_DIM
inline void Rod_shaped::find_sqr_distance_to_rod(vec2 &rr, double &dist2, double &t)
{
  double drx=rr.x-r.x, dry=rr.y-r.y ;
  t=(drx*n.x+dry*n.y) ;   
  if (t>d) t=d ;
  if (t<-d) t=-d ;  
  dist2=(SQR(drx-t*n.x)+SQR(dry-t*n.y)) ;
}


inline double sqr_distance_from_point_to_line(vec &r,  vec &p1, vec &c, double l2, double &t) 
// r=the point from where to calculate the distance
// p1 = start point of the line, p2 = end point, l2=line length
// returns the squared distance and t
{
  vec d=r-p1 ;
  t=scalar(d,c) ; t/=l2 ;  
  if (t<0 || t>1) return 1e6 ;
  d-=c*t ;
  return squared(d) ;
}

void find_forces_SS(vec &r1, vec &r2,  vec &f1);

void force_from_substrate(double m, vec2 r,vec &f)
{
  f.zero() ;
	
	const double fff=Rod_shaped::E_bact ;
	if (r.x>width/2-d0/2) f.x-=fff*(r.x-(width/2-d0/2)) ;
	if (r.x<-width/2+d0/2) f.x-=fff*(r.x-(-width/2+d0/2)) ;
	if (r.y>height/2-2*d0) f.y-=fff*(r.y-(height/2-2*d0)) ;
	
	//if (r.y>ysurface(r.x)) f.y-=fff*(r.y-ysurface(r.x)) ;
#ifdef BOTTOM_SURFACE
	vec p;
	double d2=find_sqr_dist_to_curve(r,p) ;
	if (d2<d0*d0) f-=p*fff*(sqrt(d2)-d0) ;
	else if (r.y-d0/2>ysurface(r.x)) f.y-=fff*(r.y-d0/2-ysurface(r.x)) ;
#endif
}

#else
  #error dimensionality not defined
#endif


#ifdef STICKY
// calculates the force (fx1,fy1,fz1) 
// between a pair of spheres at (rx1,ry1,rz1) and (rx2,ry2,rz2) having velocities (vx1,vy1,vz1) and (vx2,vy2,vz2)
// the force is applied to sphere 1
#ifdef FRICTION
void find_forces_SS(double rx1,double ry1,double rx2,double ry2,  double &fx1, double &fy1,  double vx1,double vy1,double vx2,double vy2)
#else
void find_forces_SS(vec &r1, vec &r2,  vec &f1)
#endif
{
  const double delta=0.05, rng=1.+2*delta, fn0=STICKY*Rod_shaped::E_bact*(1./(delta*delta*delta*delta)) ; 
  const double d20=d0*d0,drng20=d20*rng*rng ;

  vec n=r2-r1 ; 
  double n2=squared(n) ;
  #ifdef THREE_DIM
    if (n2>drng20) { f1=vec3(0,0,0) ; return ; }
	#else
  	if (n2>drng20) { f1=vec2(0,0) ; return ; }
	#endif
  double ln=sqrt(n2) ;
  n/=ln ;
  double fn ; 
  if (n2<=d20) {
    fn=Rod_shaped::pow15[int(Rod_shaped::powmax15*(d0-ln)/d0)] ; //pow(d0-n,1.5) ;
#ifdef FRICTION
		not implemented properly
    double ft=10000*pow05[int(powmax05*(d0-n)/d0)] ; //(d0-n,0.5) ;
    double dvx=vx2-vx1, dvy=vy2-vy1, sprod=(dvx*nx+dvy*ny) ;
    double vtx = dvx - sprod*nx ;
    double vty = dvy - sprod*ny ;
    fx1=-nx*fn + vtx*ft ; fy1=-ny*fn + vty*ft ;
#else
    f1=n*(-fn) ;
//    fx1=-nx*fn ; fy1=-ny*fn ;
#endif
  } else {
    //fn=-STICKY*rep_max*exp(-500*SQR(d0*1.125-n)) ; ; 
    double x=ln/(d0) ; 
    fn=-fn0*(x-1)*(x - 1)*(1 + 2*delta - x)*(1 + 2*delta - x) ;    
    //fx1=-nx*fn  ; fy1=-ny*fn  ;
    f1=n*(-fn) ;
  }
}
#else
// calculates the force (fx1,fy1,fz1) 
// between a pair of spheres at (rx1,ry1,rz1) and (rx2,ry2,rz2) having velocities (vx1,vy1,vz1) and (vx2,vy2,vz2)
// the force is applied to sphere 1
#ifdef FRICTION
void find_forces_SS(vec &r1, vec &r2,  vec &v1, vec &v2, vec &f1)
#else
void find_forces_SS(vec &r1, vec &r2,  vec &f1)
#endif
{
  const double d20=d0*d0 ;
  vec n=r2 ; n-=r1 ;
  double n2=squared(n) ;
  if (n2>d20) { f1.zero() ; return ; }
  double ln=sqrt(n2) ;
  n/=ln ;
//  nx/=n ; ny/=n ; nz/=n ;
  double fn=Rod_shaped::pow15[int(Rod_shaped::powmax15*(d0-ln)/d0)] ; //pow(d0-n,1.5) ;
#ifdef FRICTION
  double ft=fr_coeff*Rod_shaped::pow05[int(Rod_shaped::powmax05*(d0-ln)/d0)] ; //(d0-n,0.5) ;
  vec dv=v2-v1 ;
  vec vt = dv - n*scalar(dv,n) ;
  f1=n*(-fn) + vt*ft ;
//  if (fabs(norm(vt)*ft)>fabs(fn)) err(" ",ft) ;//err("x",norm(vt)*ft/fn) ;
//  fx1=-nx*fn + vtx*ft ; fy1=-ny*fn + vty*ft ;
#else
  f1=n*(-fn) ;
//  fx1=-nx*fn  ; fy1=-ny*fn  ; fz1-=nz*fn ;
#endif
/*  Lines l ;
  l.r1=r1 ; l.r2=r1+vt*ft ;
  lines.push_back(l) ;
  l.r1=r2 ; l.r2=r2-vt*ft ;
  lines.push_back(l) ;
*/
}
#endif


#ifdef ADHESION
const double max_strain=0.2 ;
const double spring_constant=1e6 ;
//const float thalf_adh=0.1 ; // time to half-adhesion
inline void force_adhesion(vec r1, vec r2, vec &df)
{
	df=(r1-r2)*spring_constant ;
  //err("x",norm(df)) ;
	//if (norm(df)>1e5) err("x",norm(df)) ;
}
#endif

//extern double ycm ;
extern float mean_front_thickness, mean_front_roughness,_fixed_roughness ;

void Rod_shaped::update_acceleration_i()
{
  vec f1,f2, df, r1=n*d,r2=n*(-d) ;

  force_from_substrate(mass,r1+r,f1) ;
  force_from_substrate(mass,r2+r,f2) ;
#ifdef ADHESION
	if (!f1.iszero() && radh1.size()==0/* && r.x<0*/) { radh1.push_back(r1+r) ; tadh1.push_back(tt) ; }
	if (!f2.iszero() && radh2.size()==0/* && r.x<0*/) { radh2.push_back(r2+r) ; tadh2.push_back(tt) ; }

	for (int i=0;i<radh1.size();i++) {
		force_adhesion(radh1[i],r1+r,df) ;
		//if (norm(df)>1e3) err("x",norm(df)) ;
		//printf("%lf ",norm(df)) ;
		double maxf=max_strain*spring_constant ; //*(tt-tadh1[i]+0.2*thalf_adh)/(tt-tadh1[i]+thalf_adh) ;
		if (norm(df)<maxf) f1+=df ;
		else {
			radh1[i]=radh1[radh1.size()-1] ; radh1.pop_back() ;
			tadh1[i]=tadh1[tadh1.size()-1] ; tadh1.pop_back() ;
			//err("x",norm(df)) ;
		}
	}
	for (int i=0;i<radh2.size();i++) {
		force_adhesion(radh2[i],r2+r,df) ;
		double maxf=max_strain*spring_constant ; //*(tt-tadh2[i]+0.2*thalf_adh)/(tt-tadh2[i]+thalf_adh) ;
		if (norm(df)<maxf) f2+=df ;
		else {
			radh2[i]=radh2[radh2.size()-1] ; radh2.pop_back() ;
			tadh2[i]=tadh2[tadh2.size()-1] ; tadh2.pop_back() ;
			//err("x") ;
		}
	}
#endif
  f+=f1 ; f+=f2 ; 

  torque+=cross(r1,f1) ;
  torque+=cross(r2,f2) ;

#ifdef TWO_DIM
  if (ext_torque>0) { double qq=phi-ext_torque_angle ; torque+=ext_torque*cos(qq-floor(qq/M_PI)*M_PI) ; }
#endif

//  f.y-=_flatten*r.y ;

//  if (_fixed_roughness>0 && mean_front_roughness>_fixed_roughness) f.y-=500*(mean_front_roughness-_fixed_roughness)*r.y ; 
  
}


void find_sqr_distance_lines(Rod_shaped *bi,Rod_shaped *bj,double &dij,double &ti, double &tj) 
{
  vec rij= bi->r, ni=bi->n, nj=bj->n, r1,r2,f1 ;
  rij-=bj->r ;
//  double rijx=bi->r.x-bj->r.x, rijy=bi->r.y-bj->r.y ; //, rij.z=bi->r.z-bj->r.z ; 
//  double nix=bi->n.x, niy=bi->n.y, /*niz=bi->n.z, */njx=bj->n.x, njy=bj->n.y ; //, njz=bj->n.z ;
  double ninj=scalar(ni,nj), w=ninj*ninj-1 ;
//  double ninj=nix*njx+niy*njy, w=ninj*ninj-1 ; 
  if (w==0) { ti=tj=0 ; }
  else {
    double snirij=scalar(ni,rij), snjrij=scalar(nj,rij) ;
//    double snirij=nix*rijx+niy*rijx, snjrij=njx*rijx+njy*rijy ;
    ti=(snirij-ninj*snjrij)/w ;
    tj=(-snjrij+ninj*snirij)/w ;
    if (ti>bi->d) ti=bi->d ;
    if (ti<-bi->d) ti=-bi->d ;
    if (tj>bj->d) tj=bj->d ;
    if (tj<-bj->d) tj=-bj->d ;
  }
//  nix*=ti ; niy*=ti ; njx*=tj ; njy*=tj ;
  ni*=ti ; nj*=tj ;
  vec d=rij ; d+=ni ; d-=nj ;
//  double dx=rijx+nix-njx, dy=rijy+niy-njy ; //, dz=rij.z+ni.z-nj.z ;
  dij=scalar(d,d) ;    
//  dij=dx*dx+dy*dy ; //+dz*dz ;
}

void Rod_shaped::add_abs_forces(vec f) 
{
//  fmag+=fabs(f.z) ;
  fmag+=manhattan(f) ; 
}

// increses accelerations of interacting rod-shaped cells i,j
//void Rod_shaped::update_accelerations_ij_rods(int i,int j) 
double Rod_shaped::update_accelerations_ij_rods(Rod_shaped *bi,Rod_shaped *bj) 
{
//  if (bi==bj) err("bi==bj") ;
//  Rod_shaped *bi=bact[i],*bj=bact[j] ;
  static vec r1,r2,f1,ri1,ri2,rj1,rj2,bin,bjn,dr1,dr2,v1,v2 ;
/*  ri1=bi->r + bi->n * bi->d ;
  ri2=bi->r - bi->n* bi->d ;
  rj1=bj->r + bj->n* bj->d ;
  rj2=bj->r - bj->n* bj->d ;
*/
  bin=bi->n ; bin*=bi->d ;
  bjn=bj->n ; bjn*=bj->d ;
  ri1=bi->r ; ri1+= bin ;
  ri2=bi->r ; ri2-= bin ;
  rj1=bj->r ; rj1+= bjn ;
  rj2=bj->r ; rj2-= bjn;

  double di1, di2, dj1, dj2, ti1, ti2, tj1, tj2, dij, ti, tj, dmin ;
//  const double d02=d0*d0 ;
  bj->find_sqr_distance_to_rod(ri1,di1,ti1) ; //if (di1<d02) goto g1 ;
  bj->find_sqr_distance_to_rod(ri2,di2,ti2) ; //if (di2<d02) goto g1 ; 
  bi->find_sqr_distance_to_rod(rj1,dj1,tj1) ; //if (dj1<d02) goto g1 ; 
  bi->find_sqr_distance_to_rod(rj2,dj2,tj2) ; //if (dj2<d02) goto g1 ; 
//  if (di1<d02 || di2<d02 || dj1<d02 || dj2<d02) { dij=10 ; goto g1 ; }
  find_sqr_distance_lines(bi,bj,dij,ti,tj) ;
//  dij=10 ;

  if (dij<1e-6) return 1e-6 ; // BW this should not happen!
//  err("dij",dij) ;
  int sgn=0 ;
  Rod_shaped *b1=NULL, *b2=NULL ;
  if (di1<=di2 && di1<=dj1 && di1<=dj2 && di1<=dij) { dmin=di1 ; b1=bi ; b2=bj ; r1=ri1 ; r2=bj->n ; r2*=ti1 ; r2+=bj->r ; sgn=1 ; goto g2 ; }
  if (di2<=di1 && di2<=dj1 && di2<=dj2 && di2<=dij) { dmin=di2 ; b1=bi ; b2=bj ; r1=ri2 ; r2=bj->n ; r2*=ti2 ; r2+=bj->r ; sgn=-1 ; goto g2 ; }
  if (dj1<=di1 && dj1<=di2 && dj1<=dj2 && dj1<=dij) { dmin=dj1 ; b1=bj ; b2=bi ; r1=rj1 ; r2=bi->n ; r2*=tj1 ; r2+=bi->r ; sgn=1 ; goto g2 ; }
  if (dj2<=di1 && dj2<=di2 && dj2<=dj1 && dj2<=dij) { dmin=dj2 ; b1=bj ; b2=bi ; r1=rj2 ; r2=bi->n ; r2*=tj2 ; r2+=bi->r ; sgn=-1 ; goto g2 ; }
  if (dij<=di1 && dij<=di2 && dij<=dj1 && dij<=dj2) { dmin=dij ; b1=bi ; b2=bj ; r1=bi->r ; r1+=bi->n*ti ; r2=bj->r ; r2+=bj->n*tj ; goto g2 ; } 

  if (b1==NULL || b2==NULL || r1==r2) {
    char txt[100] ;
    txt[0]='!' ; txt[1]=0 ;
    for (int i=0;i<bact.size();i++) if (bj==bact[i]) {
        if (isnan(bj->r.x) || isnan(bj->r.y)) err("isnan_bj") ;
        if (isnan(bi->r.x) || isnan(bi->r.y)) err("isnan_bi") ;
  	    sprintf(txt, "xyi=(%lf %lf) xyj=(%lf %lf) %lf %lf %lf %lf %lf xxx %lf %lf %lf %lf %lf %lf",
          bj->r.x,bj->r.y,bi->r.x,bi->r.y,di1,di2,dj1,dj2,dij,ti1,ti2,tj1,tj2,ti,tj);
      }
//    sprintf(txt,"%lf",dij) ;
	  err(txt) ;
  }

g2:
  dr1=r1 ; dr1-=b1->r ; dr2=r2 ; dr2-= b2->r ;
#ifdef FRICTION
  v1=b1->v ; add_cross(v1,b1->omega,dr1) ;
  v2=b2->v ; add_cross(v2,b2->omega,dr2) ;
  find_forces_SS(r1, r2, v1, v2, f1) ;    
//    vx[m] -sn[m]*vphi[m]*d[m]*sgn +cs[m]*vd[m]*sgn,  vy[m] +cs[m]*vphi[m]*d[m]*sgn +sn[m]*vd[m]*sgn, 
//    vx[n] -sn[n]*vphi[n]*t +cs[n]*vd[n]*t/d[n],  vy[n] +cs[n]*vphi[n]*t +sn[n]*vd[n]*t/d[n]) ; 
#else
  find_forces_SS(r1, r2, f1) ; 
//  if (isnan(f1.x) || isnan(f1.y)) err("ij_f1",bi->r.y-bj->r.y);
#endif
  nointeracting++ ;
  #ifdef JUNCTIONS
		if (dmin<d0 && ttt%ttt_junctions==0 && is_junction(b1,b2)==0) { // creation of new junctions every 100 steps
			 add_junction(b1,b2,r1,r2) ;
   	}
//   	for (int i=0;i<b1->junctions.size();i++) force_from_junction(  force and torque must be included ! 
  #endif


  if (f1.iszero()) return dmin;


  b1->f+=f1 ; b1->add_abs_forces(f1) ;
  b2->f-=f1 ; b2->add_abs_forces(f1) ;

//  vec cr=(r1+r2)*0.5, rb1=cr-b1->r, rb2=cr-b2->r ;
//  b1->torque+=cross(rb1,f1) ;
//  b2->torque-=cross(rb2,f1) ;
  add_cross(b1->torque,dr1,f1) ;
  subs_cross(b2->torque,dr2,f1) ;


// conjugation
  if (conjugation_rate>0 && _drand48()<dt*conjugation_rate) {
    if (b1->type==0 && (b2->type==1 || b2->type==2)) b1->type=2 ;
    if ((b1->type==1 || b1->type==2) && b2->type==0) b2->type=2 ;
  }

  return dmin ;
}

void Rod_shaped::recalc_mass_and_I() 
{
  const double q=1.*r0*r0*M_PI ;
  fmagav+=dt*(fmag-fmagav) ;
  fmag=0 ;
//    double q=1.*r0*r0*M_PI ; 
//    if (bact[i]->type==1) q*=0.1 ;
  mass=q*(2*d0*d+(4/3.)*r0*d0) ; // mass in pg (1e-15 kg)
//    double d=bact[i]->d ;
  double a=1.5*d*r0*r0+0.533*r0*r0*r0+1.333*d*d*r0+0.666*d*d*d, b=r0*r0*(d+(8./15)*r0) ;
#ifdef THREE_DIM
  mom_inertia=vec(a,a,b)*q ;  // moment of inertia, in pg*um^2
#else
  mom_inertia=a*q ;
#endif
}

void Rod_shaped::clear_forces() 
{
  f.zero() ; 
#ifdef THREE_DIM
  torque.zero() ; // initially we set all accelerations to zero
#else
  torque=0 ;
#endif
  vis=-1 ;
}

#ifdef CLOSEST
  void Rod_shaped::find_acceleration(int i)
  {
    int j,k,m;
    if (active<0) return ; // do not calculate forces for inactive cells (which do not grow and are deep in the colony)
    this->update_acceleration_i() ;
    for (m=0;m<closest.size();m++) {
  /*    for (int n=m+1;n<closest.size();n++) if (closest[m]==closest[n]) {
        err("x",int(closest.size())) ; 
        //err("x",closest[1]) ;
      }*/
      j=closest[m] ; 
      if (j>i) { // && bact[j]->vis<i) {
  //      bact[j]->vis=i ;
        double dmin=Rod_shaped::update_accelerations_ij_rods(bact[i],bact[j]) ; 
        if (dmin>SQR(d0*1.5)) { // remove neighbour if calculated minimal distance too large
          closest[m]=closest[closest.size()-1] ; closest.pop_back() ; m-- ;
  /*        int n ; 
          for (n=0;n<bact[j]->closest.size();n++) if (bact[j]->closest[n]==i) { 
              bact[j]->closest[n]=bact[j]->closest[bact[j]->closest.size()-1] ; bact[j]->closest.pop_back() ; 
              break ;
            }*/
  //        noclosest-- ;
        }
      }
    }
  
  }
#else   // if not using the proximity list then....
  void Rod_shaped::find_acceleration(int i)
  {
    int j,k;
    if (active<0) return ; // do not calculate forces for inactive cells (which do not grow and are deep in the colony)
    this->update_acceleration_i() ;

    int bxn, byn, *box ;  
    if (bi==-1 || bx==-1 || by==-1) { bxn=xtobox(r.x) ; byn=ytobox(r.y) ; } else { bxn=bx ; byn=by ; }
    for (int k=0;k<=8;k++) {
      box=&(boxes[(bxn+kx[k]+nbox)%nbox][(byn+ky[k]+nbox)%nbox][0]) ;
      int nmax=boxes[(bxn+kx[k]+nbox)%nbox[(byn+ky[k]+nbox)%nbox].size() ;
      for (int n=0;n<nmax;n++) {
        int j=box[n] ;


        if (j>i) { // && bact[j]->vis<i) {
  //      bact[j]->vis=i ;
          double dmin=Rod_shaped::update_accelerations_ij_rods(bact[i],bact[j]) ; 
        }
        
      }
  
    }
  }
#endif

void Rod_shaped::one_Newton_step(int thisi) {
}


void Rod_shaped::one_Euler_step(int thisi) {
#ifdef THREE_DIM
  if (active<1) { f.zero() ; v.zero() ; torque.zero() ; omega.zero() ; } // no acceleration
#else
  if (active<1) { f.zero() ; v.zero() ; torque=0 ; } // no acceleration
#endif

  double zetalin=zeta ; 


// anisotropic friction
  if (ANISOTROPIC_FRICTION!=0) {
    double kper=mass*zetalin*(1-ANISOTROPIC_FRICTION), kpar=mass*zetalin/(1-ANISOTROPIC_FRICTION) ;
    double denom=kpar*kper ;
    v.x=(f.x*kper*n.x*n.x-f.y*kpar*n.x*n.y+f.y*kper*n.x*n.y+f.x*kpar*n.y*n.y)/denom ;
    v.y=(f.y*kpar*n.x*n.x-f.x*kpar*n.x*n.y+f.x*kper*n.x*n.y+f.y*kper*n.y*n.y)/denom ;
  } else {
    v=f/(mass*zetalin) ; // um/h
  }
  
  vec dr=v*dt ;
  vec rrr=f ; 
  r+=dr ; 
  if (manhattan(dr)>_maxdisp) _maxdisp=manhattan(dr) ;
  if (isnan(r.x) || isnan(r.y)) { noncriterr=1 ; return ; } //err("isnan_euler") ;

  vav=vav*(1-30*dt) ;
  vav+=v*dt*30 ;

#ifdef THREE_DIM
  double a=mom_inertia.x, b=mom_inertia.z ;
  double csp2=csp*csp, snp2=snp*snp, cst2=cst*cst, snt2=snt*snt, cst4=cst2*cst2, snt4=snt2*snt2 ;
  double m11=(b*cst4*snp2 + a*csp2*snt2 + b*cst2*(csp2 + 2*snp2*snt2) +  b*snp2*snt4)/(a*b) ;
  double m12=((a - b)*snt2*csp*snp)/(a*b) ;
  double m13=((a - b)*csp*cst*snt)/(a*b) ;
  double m22=(b*csp2*cst4 + a*snp2*snt2 + b*cst2*(snp2 + 2*csp2*snt2) + b*csp2*snt4)/(a*b) ;
  double m23=((a - b)*cst*snp*snt)/(a*b) ;
  double m33=(a + b + (a - b)*(cst2-snt2))/(2*a*b) ;

  omega.x=m11*torque.x + m12*torque.y + m13*torque.z ;
  omega.y=m12*torque.x + m22*torque.y + m23*torque.z ;
  omega.z=m13*torque.x + m23*torque.y + m33*torque.z ;

  vec dn=cross(omega,n)*dt/zeta ;
  n+=dn ;
  if (manhattan(dn)>_maxdisp) _maxdisp=manhattan(dn) ;
  
  normalize(n) ;
  phi=carpol(n.x,n.y) ;
  theta=M_PI/2-carpol(sqrt(n.x*n.x+n.y*n.y),n.z) ;
              
  csp=cos(phi) ; snp=sin(phi) ;
  cst=cos(theta) ; snt=sin(theta) ;
#else // 2d case
  double dphi = dt*torque/(zeta*mom_inertia) ;
  if (d*dphi>_maxdisp) _maxdisp=d*dphi ;
  phi+=dphi ; 
  phi=fmod(phi,2*M_PI) ;
  n.x=csp=cos(phi) ; n.y=snp=sin(phi) ;
#endif

  if (_maxdisp>d0) { noncriterr=1 ; return ; } //err("jump!") ;

#ifdef RECTBOX
  if (r.x>width/2-0.25*d0) r.x=width/2-0.25*d0 ;
  if (r.x<-width/2+0.25*d0) r.x=-width/2+0.25*d0 ;
#endif

  dr=r+n*d - rold1 ; if (squared(dr)>SQR(0.3*d0)) nn_needs_update=1 ; //{ this->find_closest(thisi) ; return ; }
  dr=r-n*d - rold2 ; if (squared(dr)>SQR(0.3*d0)) nn_needs_update=1 ; //this->find_closest(thisi) ; 

}

double * conc(double x, double y)
{
  if (x<-1e6 || y<-1e6 || x>1e6 || y>1e6) err("infconc") ;
  int a=xtotab(x), b=ytotab(y) ;
  if (a<0 || a>=tabx) err("x outside tabx",a) ;
  if (b<0 || b>=taby) err("y outside taby",a) ;  
  return &(cc[a][b]) ;
}

#ifdef MONOD
double uptakefun(double c) // Monod uptake function
{
  if (c>0) return c/(c+1) ; else return 0 ;
}
#else
double uptakefun(double c) // constant uptake function (for c>0)
{
  if (c>0) return 1. ; else return 0 ;
}
#endif

void Rod_shaped::one_growth_step(int thisi) 
{
//  int a=xtotab(r.x), b=ytotab(r.y) ;
//  if (a<0 || a>=tabx || b<0 || b>=taby) err("ab") ;
//     float gr=1-SQR(1.2*b/taby) ; if (gr<0) gr=0 ;
//    float gr=growth_rate*fitness[i]*(cc[a][b]>cutoff_c?cc[a][b]:0) ; // the growth rate is proportional to the concentration of food
#ifdef NUTRIENT_LIMITED
  #ifdef MAX_ELONG_RATE
    double *c=conc(r.x,r.y) ;
    float gr=0.5*growth_rate*fitness*uptakefun((*c)) ; if ((*c)<cutoff_c) gr=0 ; 
    uptake=eat_rate*fitness ; if (uptake<0) uptake=0 ; 
  #else
    float area=(2*d0*d+2/3.*d0*d0) ;
  //  if (r.z>0.5) volume=0 ;
    double *c=conc(r.x,r.y) ;
    float gr=area*growth_rate*fitness*uptakefun((*c)) ; /*  (*c)/((*c)+1) ; */if ((*c)<cutoff_c) gr=0 ; 
    uptake=area*eat_rate*fitness ; if (uptake<0) uptake=0 ; 
  #endif
  (*c)-= uptake*uptakefun((*c))*dt/(dx0*dx0) ; if ((*c)<0) (*c)=0 ; // decrease concetration of food
//  err("x",(*c)) ;
#else
  float gr=0.5*growth_rate*fitness ;
#endif
//        if (gr>0) gr=growth_rate*fitness[i] ;
  // cc[a][b] corresponds to the box in which the cell is

/*#ifdef WASTE_LIMITED
  float gr=growth_rate*fitness[i]*(1-cc[a][b]) ; if (gr<cutoff_c) gr=0 ; 
  cc[a][b]+= waste_rate*fitness[i]*(1-cc[a][b])*dt ; if (cc[a][b]<0) cc[a][b]=0 ; // increase concetration of waste
#endif
*/


  if (active>=1) { // only active cells can grow and replicate

    vd=gr ;
    d+=dt*vd ; // increase the length of the rod
    
// BW this may still be necessary!!    if (d/dold>1.3) nn_needs_update=1 ; //this->find_closest(thisi) ; 

// BW phi[i]+=dt*vd[i] ; 
// BW x[i]+=2*cs[i]*dt*vd[i]-d[i]*dt*vd[i]*sn[i] ; y[i]+=2*sn[i]*dt*vd[i]+d[i]*dt*vd[i]*cs[i] ; // BW --- "pedalling"
#ifdef SLIME
    slime[a][b]+=dt ;  // slime is on
#endif  
    

//        if (d_final[i]>0.2 && d[i]>=d_final[i]*0.9) d[i]=d_final[i]*0.9 ; // stop growing if bigger than some critical value
    // this has an effect only if cells are allowed to be longer than 4*r0.
    
    if (tt>0 && d>=d_final) { // birth of a new cell
      // position of a new cell is calculated, the old cell is shifted
      // so that after replication the two cells touch each other with circular
      // caps as if the mother cell divided in half
      Rod_shaped *nb=new Rod_shaped(type) ;
      bact.push_back(nb) ;
      
      nb->r=r-n*(r0+d)*0.5 ;
      r+=n*(r0+d)*0.5 ;
      nb->v=v ; 
      d=(d-r0)/2.0 ; nb->d=d ; if (d<0) err("d<0 in one_growth_step()") ;
      nb->phi=phi ; 
      nb->fitness=fitness ; nb->active=1 ; nb->prev_number=number ;
      nb->pole1=pole1+1 ; pole2++ ; nb->pole2=0 ; pole1=0 ;
      nb->vd=vd ;
#ifdef MUTATION_MODE
      if (MUTATION_MODE==1 && _drand48()<mut_prob) { // mutation
#ifndef COLOR_MUTATIONS
        nb->type=1 ; 
#else
        nb->type=(++no_of_mutations) ; 

#endif
        nb->fitness=1+fit_adv ;
      }
      if (MUTATION_MODE==2) {
        if (mut_occurred==false && _drand48()<mut_prob) { // mutation
          nb->type=1 ; nb->fitness=1+fit_adv ;
          mut_occurred=true ; // only one mutation and then none
          just_mutated(bact.size()-1) ;
        }      
      }
#endif
/*        if (d_final[i]==0.2) {
        d_final[nn]=0.2 ;
        if (drand48()<p_filam) d_final[nn]=0.1*max_len ;
      } else {
        d_final[i]=d_final[nn]=0.2 ;
      }*/

      double qq=0.002 ;      
      phi+=(_drand48()-0.5)*qq ; nb->phi+=(_drand48()-0.5)*qq ; // some randomness introduced

#ifdef THREE_DIM
      nb->theta=theta ;       
      theta+=(_drand48()-0.5)*qq ; nb->theta+=(_drand48()-0.5)*qq ; // some randomness introduced
#endif

#ifdef RANDOMIZE_DIRS
      phi=_drand48()*M_PI*2 ; nb->phi=_drand48()*M_PI*2 ; // randomized direction
#ifdef THREE_DIM
      theta=_drand48()*M_PI ; nb->theta=_drand48()*M_PI ; // randomized direction
#endif
#endif
      this->recalc_cspt() ;
      nb->recalc_cspt() ;
      
      this->nn_needs_update=1 ;
      nb->nn_needs_update=1 ;
      
//      double q=drand48() ; q*=q ;//*q ;
//      d_final=d_final*(0.9+drand48()*0.2) ;
//      nb->d_final=d_final*(0.9+drand48()*0.2) ;
    }
  }
  
}


void Rod_shaped::recalc_cspt() 
{
#ifdef THREE_DIM
  csp=cos(phi) ; snp=sin(phi) ; cst=cos(theta) ; snt=sin(theta) ; n=vec(csp*snt,snp*snt,cst) ; 
#else
  csp=cos(phi) ; snp=sin(phi) ; n=vec(csp,snp) ; 
#endif
}

#ifdef CLOSEST
  #ifdef NON_ACTIVE_STATIC
int Rod_shaped::change_state() 
{
  int ret=0 ;
  if (active==1 && vd==0) {
    int ac=0 ;
    for (int n=0;n<closest.size();n++) if (bact[closest[n]]->vd>0) { ac=1 ; break ; }
    if (ac==0) { active=0 ; ret=1 ; }
  }
  if (active==0) {
    int ac=0 ;
    for (int n=0;n<closest.size();n++) if (bact[closest[n]]->active>0) { ac=1 ; break ; }
    //int j=int(r.x+width/2) ;
    if (ac==0 && min_front_thickness>5*d0) { active=-1 ; ret=1 ; }
  }
  return ret ;
}
  #else
int Rod_shaped::change_state() { return 0 ; }   // no change to the state if !NON_ACTIVE_STATIC
  #endif
#endif
