// v1.3.2
//---------------------parameters-----------------------------

//#define CEREVISIAE
#define ECOLI 
//#define SPHERES
//#define SHORT_RAND_RODS

//#define NON_ACTIVE_STATIC // if defined, non-active cells do not move and do not interact. Works only in the expanding colony model

#define OVERDAMPED // if defined, dynamics of bacteria is overdamped, otherwise it is Newtonian

//#define MONOD // if defined, use the Monod uptake function, otherwise use constant uptake function

#define RECTBOX  // bacteria growing in a rectangular compartment of constant width, with walls on all but the top side

//#define NUTRIENT // if defined, simulate the diffusion and uptake of the nutrient

#define INITIAL_CONDITION 1 // not used for RECTBOX
// 1 = single cell
// 2 = two cells separated by some distance (set in the reset procedure)
// 3 = many cells on a ring
// 4 = many cells in a circle

//#define OBSTACLES // if defined, diffusion/bacteria can go only to certain parts of the tube

//#define CORR_DISORDER // introduces correlated disorder into the friction coefficient - not implemented yet!
#define MUTATION_MODE 5 // if defined, mutations occurs during growth. their fitness advantage is fit_adv. 
/*  modes:
  1 = cells mutate upon replication with mut. probability mut_prob
  2 = only a single cell mutates upon division
  3 = a single cell is selected at random from the first line of cells
  4 = at t=0 all cells are assigned different types from 0 to no_of_cells-1
  5 = cells are assigned two types at random
  6 = a single random cell is selected from all cells
*/

//#define REASSIGN_TYPES_N 30 // if defined, reassign diffrent types from 0 to no_of_cells-1 when the size exceeds the given value
//#define REASSIGN_TYPES_T 3 // if defined, reassign diffrent types from 0 to no_of_cells-1 after a given time

//#define COLOR_MUTATIONS // if defined, different mutations are ascribed different, random colours

#define TWO_DIM // if defined the colony is simulated in 2d which is faster than the fully 3d simulation
//#define THREE_DIM // if defined the colony is simulated as fully 3d

//#define MARK_NN_OF_0 // if defined it marks neighbours of cell 0 with different colours
//#define MARK_NONACTIVE // if defined, non-active cells are drawn as thin lines instead of rods
//#define WHITE_LINES_INSIDE // draw thin white lines inside rod-shaped cells
//#define HALF_THICKNESS // draw cells with half of the normal thickness
#define AUTOSCALE // scale so that whole colony is always visible. Ignored for TUBE

// !!!do not use when recording videos!!!
//#define TRY_INCREASE_DT // if defined, try to increase dt by 2x from time to time to speed the simulation. Does not cause any problems

//#define VERBOSE // display many more parameters on the screen in the graphical mode
//#define SAVE_FRAMES_TO_DAT // if defined, save frames to independent *.dat files (works only with main_graphics.cpp)

#define BOTTOM_SURFACE // use this to implement a non-trivial bottom surface that may be curved
#define ADHESION // adhesion between bacteria and the substrate                
//#define JUNCTIONS // create physical junctions between neighbouring cells

#ifdef __MAIN

// ------------ parameters of the simulation ----------------------

#ifdef BOTTOM_SURFACE
double bottom_period=10, bottom_ampl=0 ;
#endif

bool backup=true ; // store old configuration in case the simulation crashes and dt needs to be decreased
bool sectors=false ; // write sectors - works only for TUBE and RECTBOX
int save_cells=1+2+16+64 ; // if non-zero, save all and save positions etc. of all cells. 
                    //1==save_all | 2==save_data | 4==save_profile | 8==save_envelope
                    //16==save numbers of each type of cells
                    //32==save positions of all cells every time frame
                    //64==use average velocities when saving data
int width=100, height=110 ; //110 ; // size [um] of the simulation box, otherwise not used
double depth=3 ; // used only for 3d simulation
// width must be a multiple of boxwidth=10*d0

#ifdef THREE_DIM
double dtmax=1./(1<<15) ; // max. integration step, [h]
#else
double dtmax=1./(1<<15) ; // max. integration step, [h] (1<<16 for 80x100)
#endif
double diff_rate=20000 ; //2000. ; // rate of diffusion of nutriens [um^2/h] - 50 for most simulations
double c0=1 ; // amount of the nutrient in one box 1um x 1um
double waste_rate=0*2.6 ; // not currently implemented
float mut_prob=0*1e-2 ; // probability that a cell mutates during replication. Used only if MUTATION_MODE defined
float conjugation_rate=0.0 ; // conjugation rate (used only if non-zero)
float fit_adv=0.5 ;

#ifdef ECOLI
// E. coli
  double ANISOTROPIC_FRICTION=0 ; // if !=0, 0<a<1 means it is easier to roll. a<0 means it is easier to slide
  double r0=0.5 ; // radius of the rod [um] (bacteria are rod-shaped with 
              // semi-circular caps on both ends, each having also radius r0 )            
  double Rod_shaped::eat_rate=1.8 ; //overriden from command line // rate at which bacteria eat nutrients, [1/h] --- 3 for decay to 20%
  double Rod_shaped::cutoff_c=0.01 ; // minimal nutrient concentration below which growth ceases (previously always 0.2)
  double Rod_shaped::E_bact=1e+7 ; //4e+6 is fine for dt13, and 20um growing layer ; //3.75e+5 ; // elastic modulus of the bacterium, Pa
// BW  double Rod_shaped::growth_rate = 1.2 ; //1.65 ; // rate of growth of the wild type for cc=infty, in [um/h*pg]  
  double Rod_shaped::growth_rate = 1. ; // BW: elongation rate in [um/h]   was 4 for 37 deg C
  double Rod_shaped::max_length=1.0 ; // dimensionless. 1/2 x maximal length of the rod (without caps) before division - 2==E. Coli, 0.5=circles
  double Rod_shaped::grav=0.01 ;
  double Rod_shaped::zeta=500 ; // dynamic friction, Pa*h. 500 for fast simulations in which the density however increases towards the centre. 
                                // 10 for slow simulation with little overlap between cells
#endif

#ifdef SHORT_RAND_RODS
// short, randomized rods
  double ANISOTROPIC_FRICTION=0 ; // if !=0, 0<a<1 means it is easier to roll. a<0 means it is easier to slide
  double r0=2.0 ; // radius of the rod [um] (bacteria are rod-shaped with 
              // semi-circular caps on both ends, each having also radius r0 )            
  double Rod_shaped::eat_rate=2000; //1.8 ; //overriden from command line // rate at which bacteria eat nutrients, [1/h] --- 3 for decay to 20%
  double Rod_shaped::cutoff_c=20 ; // minimal nutrient concentration below which growth ceases
  double Rod_shaped::E_bact=1e+7 ; //3.75e+5 ; // elastic modulus of the bacterium, Pa
  double Rod_shaped::growth_rate = 2. ; // elongation rate in [um/h]  
  double Rod_shaped::max_length=0.5 ; // dimensionless. 1/2 x maximal length[um] of the rod (without caps) before division - 2==E. Coli, 0.5=circles
  double Rod_shaped::grav=1 ;
  double Rod_shaped::zeta=50 ; // dynamic friction, Pa*h. 500 for fast simulations in which the density however increases towards the centre. 
                                // 10 for slow simulation with little overlap between cells
#endif

// -------------------------------------------------------------
float _flatten=0 ; // if non-zero, a force flattening the leading edge is introduced
float ext_torque=0 ;	// if non-zero then apply external torque aligning cells
float ext_torque_angle=0*0.5*M_PI ;	// if 0 then try to align cells vertically. if M_PI/2 then horizontally
float _fixed_roughness=0 ; // if non-zero, the front is forced to never grow above this roughness, this is done by applying forces as for _flatten>0
                            // but conditional on the front becoming too rough. This effectively fixes the roughness of the front

float fmag_transition=0*5e4 ; // if non-zero, defines the force [in pN] at which the bacterium vanishes from the 2d colony (as if going to the 2nd layer)


// this is used only for WINDOWS and graphics
int winx=1000,winy=1000 ; // bitmap size. Don't change
int _stop=0 ; // if ==1, the simulation is stopped but still displayed
double _maxdisp ; // maximal displacement in one step, used to monitor convergence
int shift_cm_y=0 ; //-20 ;//  position of the cm in y-direction -- 9

#endif

//------------------------------------------------------------

#define CLOSEST // use proximity list (n.n.) to calculate interactions if agarose not modelled explicitly

