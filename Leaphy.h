/*	******************************************************************
		SWitch header file
	****************************************************************** */

/* LeaPHY and SNAP TODO list:

  1) Ensure TABU search is working correctly - very probably
  2) Debug error with some stepwise tree rearrangements (unclear what error is)
  3) Code so that multiple rearrangements can occur in one round of SNAP
  4) Exit criteria seem to not be working properly (OptObs)

  LeaPHY and SNAP ideas list:

  1) BIG IDEA: When importance sampling in tree space using SA, why not rearrange the smaller tree?
	This will allow bigger changes in topology and hopefully improve the estimation procedure.
	The natural extension of this is to SNAP every time a sequence is added, but this is likely to be too slow
*/

#ifndef __SWitch_HEADER
#define __SWitch_HEADER

#ifdef _WIN32
#pragma warning (disable:4786)
#endif

#include "tools.h"

class CData;
class CTree;

/////////////////////////////////////////////////////////////////////////////////
// General flags
/////////////////////////////////////////////////////////////////////////////////

// Developer's flag
#define DEVELOPER_VERSION 1					// Whether the current build is the developer version

// Other general flags for programs other than Leaphy
#define IS_LEAPHY 1							// Whether the program is Leaphy
#define IS_MODELCHECKER 1					// Whether the program is ModelChecker
#define IS_LNL_CALC 1						// Whether libraries are being used to do likelihood calculations

#if IS_LEAPHY == 1
#define PROG_NAME "Leaphy"
#define VERSION_NUM "1.01beta"
#elif IS_MODELCHECKER == 1
#define PROG_NAME "ModelChecker"
#define VERSION_NUM "0.1alpha"
#endif

////////////////////////////////////////////////////////////////////////
EModel string_to_model(string s);				// Transforms a string to a model number

/////////////////////////////////////////////////////////////////////////
// Program defaults
#define DEFAULT_DNA_MODEL REV
#define DEFAULT_AA_MODEL WAG
#define DEFAULT_GAMMA true
#define DEFAULT_GARBAGE false
#define DEFAULT_GAMMA_CAT 4
#define DEFAULT_INV false

// The type of stepwise addition performed
#define PROB_DO_MP_SA 0.25		// Probability of doing straight Maximum parsimony stepwise addition
#define DO_MP_HILLCLIMB	0.25		// Whether to do a maximum parsimony hill climb after MP resampling
#define DO_MP_SA 0				// Whether to do Maximum parsimony starting tree

// Whether to do tabuing
#define DO_TABU 1
#define DEFAULT_TABU_RADIUS 3	// Tabu radius of each tree when resampling from tree-space
#define ENCOURAGE_SNAP_AWAY 0	// Push away trees during SNAP -- When multiple improvements, choose the one that maximises MinRFDist to tabu treesassert(false)
#define COMP_SNAP_AWAY_2_OPT_OBS 0 // Nunber of attempts toward end (if ENCOURAGE_SNAP_AWAY = =1) that force the SNAP away
// General tree varying routine variables
enum LazyType { fullopt, lazy, randomtree };		// The type of optimisation done when performing actions
#define ALLOW_TABU_PART_OPT 1		// Allows lazy optimisation to make trees tabu...
#define LAZY_SUBTREE lazy				// Do lazy subtree optimisation
#define LAZY_SA lazy					// Do lazy stepwise addition

// Probability approximations for sitewise likelihoods
#define ALLOW_SCALE 1			// Whether likelihoods are scaled
#define P_SCALE_VAL (0.1)	// Value at which scaling takes place
#define ALLOW_FAST_CALC 0		// Whether fast likelihood computations are allowed

////////////// Sequence trimming and rearranging processes ///////////////
#define DO_SNAP 0				// Whether to perform only SNAP
#define DO_SPR_FULL 1			// Whether to perform only SPR
#define DO_IMPSAMPLE 1			// Whether to perform sampling
#define OUTPUT_TABULIST 1		// Whether to output the tabulist
#define DO_SPR_SECOND 1				// Whether to perform SPR after SNAP
#define SPR_LOWBOUND 1			// Lower bound of SPR rearrangement distance
#define SPR_UPBOUND 5			// Upper bound of SPR rearrangement distance
#define DESP_THRESHOLD 2		// Number of PSA attempts using ML before doing random rearrangements
#define PARS_DESP_THRESHOLD 5	// Number of PSA attempts using MP before doing random rearrangements
// The approximations used in the software
#define QUICK_NUMERICAL_DERIVATE 1				// [0] Takes 3 points for derivative; [1] Takes 2 points for derivative in CBaseModel::GetNumDerivative function
#define ALLOW_ANALYTICAL_BRANCH_DERIVATIVES 1  	// Whether analytic derivatives for branches are allowed
#define THOROUGH_LINE_SEARCH 5	// Number of thorough line searches to perform during beginning of optimisation
#define NUM_FAST_BRA_OPT 5		// Number of fast branch optimisations performed
#define FULL_LIK_ACC 1.0e-7		// Accuracy of full likelihood routines
#define TREESEARCH_LIK_ACC 1.0e-2	// Accuracy of tree-search likelihood
#define FULL_GTOL 1.0e-4		// Accuracy of gradient
#define PART_GTOL 1.0e-2		// Accuracy of the gradient for approximate routines
#define PART_LIK_ACC 1.0e-1		// Accuracy of approximate likelihood routines
#define RMSD_ACC 1.0e-2			// Accuracy of RMSD optimisation
#define DX 5.0E-5				// Delta(x) used for optimisation
#define PROB_DX 1.0E-4			// Delta(x) for probabilities
#define BOUND_GRAD 0			// Whether gradients should be bounded at GRAD_LIM
#define GRAD_LIM 1000			// Maximum value allowed in derivative functions
#define RMSD_GRAD_LIM 10		// Maximum value allowed in RMSD derivative functions
#define MAX_BRANCH 10			// Maximum branch length
#define MIN_OMEGA 0.05			// Minimum value of omega in codon models (required to avoid interaction between small omega and long branch lengths)
#define BIG_STEP_MAX 10			// MAX STEP when fewer than 10 iterations
#define SMALL_STEP_MAX 5		// MAX_STEP when >= 10 iterations
#define DX_CHECK 0.0001			// Left and right check for optima when gradients are playing up
#define BIG_LNL_DIFF 5			// A large difference in likelihood that is unlikely to be overcome by small optimisations
#define GOLDEN_NUMBER 0.38197		// Number for golden section search

// Parameter initialisations
#define ALLOW_KAPPA_GUESS 1		// Allows the process to guess the kappa value from the data (otherwise INITIAL_KAPPA is used)
#define INITIAL_KAPPA 2.5		// The value kappa is set to by default in models
#define INITIAL_OMEGA 0.4		// The value omega in codon models is set to by default

//////////////////////////////////////////////////////////////////////////////
// SPR approximations
// ---
#define LAZY_SPR_RATIO 1.0075
#define FULL_SPR_RATIO 1.015

//////////////////////////////////////////////////////////////////////////////
//////////////////////////// SNAP approximations //////////////////////////////
// Approximations made when assessing 4,5, and 6 species trees
#define DO_QUICK_SNAP 1			// Does quick SNAP where only one change done at a time...
#define MAX_KNOT_UNTANGLE 3		// Number of times a change can be made without improving likelihood
#define ALLOW_SNAP_SKIP 1		// Whether to allow SNAP routines to skip previously unimproved nodes
#define DELTA_STEP_SUB 0.1		// Step used to assess likelihoods in SUBTREE routines
#define DO_LEAF_SNAP 0			// Don't allow SNAP to skip nodes where 2 links are leaves
#define STEPS_BETWEEN_PAROPT 5	// Number of steps between parameter optimisation in Lazy optimisation
// NUM_ACC
#define LOOSE_OPT_NUM_ACC 4		// Number of accurate optimisations made in DoSubTree when doing loose calcs
#define FULL_OPT_NUM_ACC 4		// Number of accurate optimisations made in DoSubTree when doing accurate calcs
// PERCENT_CHANGE w.r.t. RMSD
#define LOOSE_PERCENT_CHANGE 1.05	// minimum difference in RMSD between current best and other trees
#define FULL_PERCENT_CHANGE 1.1
// Methods for choosing subsets of SNAP trees to examine
// RMSD_SUBSET -- minimum number top RMSD trees stored
#define LOOSE_RMSD_SUBSET 106
#define FULL_RMSD_SUBSET 106
// PARS_SUBSET -- minimum number top PARS trees stored
#define LOOSE_PARS_SUBSET 10
#define FULL_PARS_SUBSET 25
#define PROP_SITES_PARS 0.75		// The maximum MP score per site relative to the number of species in the tree.


/////////// Step-wise addition methods /////////////////////////////////////////
// Approximations made when performing large scale rearrangements
///////////////// Some definitions for random plucking (i) ////////////////////////
#define DEFAULT_PROB_RAN_SEQ_REM 0.33			// Probability of removing each sequence
#define PROB_DO_IMPSNAP 0.25					// Whether to do another SNAP when importance sampling
// Some definitions for random snaps
#define PROB_NODE_SNAP 0.2	// Number of randoms SNAPs performed = 1 + #internal nodes * PROB_NODE_CHANGE
#define MIN_NUM_RANSNAPS 6	// Minimum number of random SNAPs on a tree
// Degree of optimisation for SA between adding sequences
#define SHORT_OPT 5			// Max number of iterations
// Number of triplets to do full optimisation on
#define SA_OPT_NUM_ACC 4	// TEMPORARY ___ _TOREMOVE
#define SA_STORE_PER_ROUND_MP 3				// Number of candidate trees stored each round when doing SA with MP
#define SA_STORE_PER_ROUND_ML 1				// Number of candidate trees stored each round when doing SA with ML
// Some definitions for randomised stepwise addition
#define DESPERATE_SUBOPT_PROB 0.05	// When STORE_PER_ROUND == 1 and the routine becomes desperate this is the value of
									// progressively worse trees in the AddSequences routine being accepted
#define SA_SCORE 10000			// Defines the truncation of the geometric distribution
#define SCORE_FACTOR 0.005		// Defines the number that the score is *= between trees
#define MULOPT_FACTOR 0.05		// Defines the number that the score is *= when multiple optima occur.
// Definitions specifying how much optimimsation goes on between adding sequences
#define NUM_OPT_SEQ_ADD 1		// Number of rounds of MulD optimisation that occur after adding a sequence
#define DEFAULT_OPTNUM 20 
	// Number of values < DBL_EPSILON in loglikelihood function bails and returns -BIG_NUMBER

// Definition of different THMM equilibrium classes
enum ETHMM_EQM_TYPE {equ,obs,complex,complex_GC};
// Definition of different THMM hidden state classes
enum EHIDDEN_TYPE { H_none, H_same, H_diff }; // H_none = no hidden transitions (H=0); H_same = 1 hidden transition; H_diff = All seperate hidden transition
#endif
