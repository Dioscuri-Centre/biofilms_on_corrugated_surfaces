#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;
#include "functions.h"
#include "vecs.h"
#include "classes.h"
#include "params.h"


double Junction::max_strain=1.25 ; // 1.02 in Paul's simulations
double Junction::spring_constant=1e6 ;

extern vector <Junction*> junctions ;


void add_junction(Rod_shaped *b1, Rod_shaped *b2, vec p1, vec p2) 
{
  vec dr1=p1-b1->r, dr2=p2-b2->r, dp=p2-p1 ;
  double t1=scalar(dr1,b1->n)/b1->d ;
  double t2=scalar(dr2,b2->n)/b2->d ;
  //err("t1",t2) ;
  double l0=norm(dp) ;
  if (l0>1) err("lo",l0) ;
  Junction *newj=new Junction(b1,b2,t1,t2,l0) ;
  junctions.push_back(newj) ;
}

void add_junction(Rod_shaped *b, vec pb, vec psurf)  // for junctions with agar
{
  vec dr=pb-b->r, dp=pb-psurf ;
  double t=scalar(dr,b->n)/b->d ;
  double l0=norm(dp) ;
  if (l0>1) err("lo",l0) ;
  //err("t",t) ;
  Junction *newj=new Junction(b,t,psurf,l0) ;
  junctions.push_back(newj) ;
}

int is_junction(Rod_shaped *b1,Rod_shaped *b2) // very inefficient!
{
  for (int i=0;i<junctions.size();i++)
    if ((junctions[i]->b1==b1 && junctions[i]->b2==b2) || (junctions[i]->b1==b2 && junctions[i]->b2==b1)) return 1 ;
  return 0 ;
}

void Junction::add_force_from_junction() 
{
	if (b2!=NULL) { // bacteria-bacteria
		vec dr1=b1->n*t1*b1->d ;
		vec dr2=b2->n*t2*b2->d ;
		vec dr=(b2->r+dr2)-(b1->r+dr1) ;
	  double ll=norm(dr) ;
	  double ff=(ll-l0)*spring_constant ;
	  vec f1=dr*ff/ll ;
	  //if (ll/l0>max_strain) err("max strain",ll/l0) ;
	  if (ll/l0>max_strain) { // label for breaking
	    b1=NULL ;	return ;
	  }
	  
	  b1->f+=f1 ; b1->add_abs_forces(f1) ;
	  b2->f-=f1 ; b2->add_abs_forces(f1) ;
	
	  add_cross(b1->torque,dr1,f1) ;
	  subs_cross(b2->torque,dr2,f1) ;
	} else {  // with agarose
		vec dr1=b1->n*t1*b1->d ;
		vec dr=(rsurf)-(b1->r+dr1) ;
		const double freq=10 ;
	  dr.x+=0.02*rsurf.x*(cos(2*M_PI*freq*tt)-cos(2*M_PI*freq*t_creation)) ; // BW surface oscillations
	  double ll=norm(dr) ;
	  double ff=(ll-l0)*spring_constant ;
	  vec f1=dr*ff/ll ;
	  //err("x",ll) ;
	  //if (ll/l0>max_strain) err("max strain",ll/l0) ;
	  if (ll/l0>max_strain) { // label for breaking
	    b1=NULL ;	return ;
	  }
	  
	  b1->f+=f1 ; b1->add_abs_forces(f1) ;
	  add_cross(b1->torque,dr1,f1) ;
	
	}
 // err("x") ;

}


void simulate_junctions() 
{
//	if (ttt==19) err("ttt") ;
	//if (junctions.size()>100) err("size") ;
  for (int i=0;i<junctions.size();i++) {
		junctions[i]->add_force_from_junction() ;
		if (junctions[i]->b1==NULL) { // break junction
		  junctions[i]=junctions[junctions.size()-1] ; junctions.pop_back() ;
		  i-- ;
		}
  }
}


