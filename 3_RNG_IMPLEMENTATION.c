//	************************************************************************************************************************** //
//	****                                                                                                                  **** //
//	****                1st GEN. LUESCHERS LUXURY RANDOM NUMBER GENERATOR IMPLEMENTATION                                  **** //
//	****                                                                                                                  **** //
//	************************************************************************************************************************** //


struct RANLUX_STATE{
	
    unsigned long int RLTable[24];

    unsigned int cnt1; //i
    unsigned int cnt2; //j
    unsigned int step; //n
    unsigned int jump; //skip
    unsigned int mem;  //carry

    unsigned int EXTSEED;

};

static const unsigned long int CONST = 16777216;        /* 2^24 */	//two24
static const unsigned long int HEX1 = ~0x00ffffffUL;			//mask_hi
static const unsigned long int HEX2 = 0x00ffffffUL;  /* 2^24 - 1 */ 	//mask_lo

static void RANLUX_INITIALIZE (void *CURSTATE, unsigned long int EXTSEED, unsigned int LUXURY)
{
	RANLUX_STATE *MYSTATE = (RANLUX_STATE *) CURSTATE;

	MYSTATE->EXTSEED=EXTSEED;
	
	long int SEED;
	int cnt;
	
	if (EXTSEED==0) EXTSEED = 314159265;
	
	SEED=EXTSEED;
	//------------------------------------------------------------
	//	USING THE SPECIFIC IMPLEMENTATION OF F. JAMES
	//------------------------------------------------------------
	
	for( cnt=0; cnt<24; cnt++){

		unsigned long int tmp = SEED / 53668;
		SEED = 40014 * (SEED-tmp*53668)-tmp*12211;
		
		if(SEED<0) SEED+=2147483563;
			
		MYSTATE->RLTable[cnt] = SEED % CONST;

	}
	
	
	MYSTATE->cnt1=23;
	MYSTATE->cnt2=9;
	MYSTATE->step=0;
	MYSTATE->jump= LUXURY-24;
	
	if(MYSTATE->RLTable[23] & HEX1) {
		
		MYSTATE->mem=1;
		
	}else{
		
		MYSTATE->mem=0;
		
	}
	
}

static void RANLUX_INIT223 (void *CURSTATE, unsigned long int EXTSEED)
{

	RANLUX_INITIALIZE(CURSTATE,EXTSEED,223);
	
}

static void RANLUX_INIT389 (void *CURSTATE, unsigned long int EXTSEED)
{
	RANLUX_INITIALIZE(CURSTATE,EXTSEED,389);
	
}

static inline unsigned long int RANLUX_UPDATE ( RANLUX_STATE * CURSTATE)
{
	unsigned int cnt1 = CURSTATE->cnt1;
	unsigned int cnt2 = CURSTATE->cnt2;
	
	long int diff = CURSTATE->RLTable[cnt2]-CURSTATE->RLTable[cnt1]-CURSTATE->mem;
	
	if(diff&HEX1){
		CURSTATE->mem=1;
		diff &= HEX2;
	}else{
		CURSTATE->mem=0;
	}
	
	CURSTATE->RLTable[cnt1] = diff;
	
	if(cnt1==0){ cnt1=23; } else { cnt1--; }
	
	CURSTATE->cnt1=cnt1;
	
	if(cnt2==0){ cnt2=23; } else { cnt2--; }
	
	CURSTATE->cnt2=cnt2;

	return diff;
}

static inline unsigned long int RANLUX_GEN_INT (void *CURSTATE)
{
	RANLUX_STATE *MYSTATE = (RANLUX_STATE *) CURSTATE;
	
	const unsigned int jump = MYSTATE->jump;
	unsigned long int VAL = RANLUX_UPDATE(MYSTATE);

	MYSTATE->step++;

	if (MYSTATE->step == 24) {
		
		unsigned int i;
	
		MYSTATE->step = 0;
		
		for (i = 0; i < jump; i++) RANLUX_UPDATE (MYSTATE);
		
	}

	return VAL;
}


static double RANLUX_GEN_DOUBLE (void *CURSTATE)
{
  return RANLUX_GEN_INT (CURSTATE) / 16777216.0;
}


//	************************************************************************************************************************** //
//	****                                                                                                                  **** //
//	****                2nd GEN. LUESCHERS LUXURY RANDOM NUMBER GENERATOR IMPLEMENTATION FOR DOUBLE PRECISION             **** //
//	****                                                                                                                  **** //
//	************************************************************************************************************************** //

struct RANLUXD_STATE{
	
    double RLDTable[12];
    double mem;

    unsigned int cnt1; //ir
    unsigned int cnt2; //jr
    unsigned int cnt1_old; //ir_old
    unsigned int jump; //pr

    unsigned int EXTSEED;

};

static const int COUNTER_CYCLE[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0};

static const double MIN_ELEM = 1.0 / 281474976710656.0;  /* 1/2^48 */

#define RANLUX_STEP(x1,x2,i1,i2,i3)      \
          x1=RLDTable[i1] - RLDTable[i2];\
          if (x2 < 0)                    \
          {                              \
            x1-=MIN_ELEM;                \
            x2+=1;                       \
          }                              \
         RLDTable[i3]=x2

static void RANLUXD_INITIALIZE (void *CURSTATE, unsigned long int EXTSEED, unsigned int LUXURY)
{
  RANLUXD_STATE *MYSTATE = (RANLUXD_STATE *) CURSTATE;

  MYSTATE->EXTSEED=EXTSEED;
  
  int ibit, jbit, i, k, l, xbit[31];
  double x, y;

  long int SEED;

  if (EXTSEED == 0)
    EXTSEED = 1;                      /* default seed is 1 */

  SEED = EXTSEED;

  i = SEED & 0xFFFFFFFFUL;

  for (k = 0; k < 31; ++k)
    {
      xbit[k] = i % 2;
      i /= 2;
    }

  ibit = 0;
  jbit = 18;

  for (k = 0; k < 12; ++k)
    {
      x = 0;

      for (l = 1; l <= 48; ++l)
        {
          y = (double) ((xbit[ibit] + 1) % 2);
          x += x + y;
          xbit[ibit] = (xbit[ibit] + xbit[jbit]) % 2;
          ibit = (ibit + 1) % 31;
          jbit = (jbit + 1) % 31;
        }
      MYSTATE->RLDTable[k] = MIN_ELEM * x;
    }

  MYSTATE->mem = 0;
  MYSTATE->cnt1 = 11;
  MYSTATE->cnt2 = 7;
  MYSTATE->cnt1_old = 0;
  MYSTATE->jump = LUXURY;
}

static void RANLUXD_INIT202 (void *CURSTATE, unsigned long int EXTSEED)
{
	
	RANLUXD_INITIALIZE(CURSTATE,EXTSEED,202);
	
}

static void RANLUXD_INIT397 (void *CURSTATE, unsigned long int EXTSEED)
{
	RANLUXD_INITIALIZE(CURSTATE,EXTSEED,397);
	
}

static inline void RANLUXD_UPDATE (RANLUXD_STATE *CURSTATE)
{
	
  int k, kmax;
  double y1, y2, y3;

  double *RLDTable = CURSTATE->RLDTable;
  double mem = CURSTATE->mem;
  unsigned int cnt1 = CURSTATE->cnt1;
  unsigned int cnt2 = CURSTATE->cnt2;

  for (k = 0; cnt1 > 0; ++k)
    {
      y1 = RLDTable[cnt2] - RLDTable[cnt1];
      y2 = y1 - mem;
      if (y2 < 0)
        {
          mem = MIN_ELEM;
          y2 += 1;
        }
      else
        {
          mem = 0;
        }
      RLDTable[cnt1] = y2;
      cnt1 = COUNTER_CYCLE[cnt1];
      cnt2 = COUNTER_CYCLE[cnt2];
    }

  kmax = CURSTATE->jump - 12;

  for (; k <= kmax; k += 12)
    {
      y1 = RLDTable[7] - RLDTable[0];
      y1 -= mem;

      RANLUX_STEP (y2, y1, 8, 1, 0);
      RANLUX_STEP (y3, y2, 9, 2, 1);
      RANLUX_STEP (y1, y3, 10, 3, 2);
      RANLUX_STEP (y2, y1, 11, 4, 3);
      RANLUX_STEP (y3, y2, 0, 5, 4);
      RANLUX_STEP (y1, y3, 1, 6, 5);
      RANLUX_STEP (y2, y1, 2, 7, 6);
      RANLUX_STEP (y3, y2, 3, 8, 7);
      RANLUX_STEP (y1, y3, 4, 9, 8);
      RANLUX_STEP (y2, y1, 5, 10, 9);
      RANLUX_STEP (y3, y2, 6, 11, 10);

      if (y3 < 0)
        {
          mem = MIN_ELEM;
          y3 += 1;
        }
      else
        {
          mem = 0;
        }
      RLDTable[11] = y3;
    }

  kmax = CURSTATE->jump;

  for (; k < kmax; ++k)
    {
      y1 = RLDTable[cnt2] - RLDTable[cnt1];
      y2 = y1 - mem;
      if (y2 < 0)
        {
          mem = MIN_ELEM;
          y2 += 1;
        }
      else
        {
          mem = 0;
        }
      RLDTable[cnt1] = y2;
      cnt1 = COUNTER_CYCLE[cnt1];
      cnt2 = COUNTER_CYCLE[cnt2];
    }
    
  CURSTATE->cnt1 = cnt1;
  CURSTATE->cnt1_old = cnt1;
  CURSTATE->cnt2 = cnt2;
  CURSTATE->mem = mem;
  
}


static double RANLUXD_GEN_DOUBLE (void *CURSTATE)
{
	RANLUXD_STATE *MYSTATE = (RANLUXD_STATE *) CURSTATE;
	
	int curcnt = MYSTATE->cnt1;

	MYSTATE->cnt1 = COUNTER_CYCLE[curcnt];

	if (MYSTATE->cnt1 == MYSTATE->cnt1_old)
	RANLUXD_UPDATE (MYSTATE);

	return MYSTATE->RLDTable[MYSTATE->cnt1];	
	
}


  
static inline unsigned long int RANLUXD_GEN_INT (void *CURSTATE)
{
  return RANLUXD_GEN_DOUBLE (CURSTATE) * 4294967296.0;
}

//	************************************************************************************************************************** //
//	****                                                                                                                  **** //
//	****                   THE INFRASTRUCTURE TO ACCOMODATE SEVERAL DIFFERENT RANDOM NUMBER GENERATORS (LIKE GSL)         **** //
//	****                                                                                                                  **** //
//	************************************************************************************************************************** //


struct TLP_ALGORITHM {
	
    int ALGORITHM;
    unsigned long int MAX_INT;
    unsigned long int MIN_INT;
    int  STATE_SIZE;
    
    void (*SET_RND) (void *CURSTATE, unsigned long int EXTSEED);
    unsigned long int (*GET_INT) (void *CURTSTATE);
    double (*GET_DOUBLE) (void *CURSTATE);
	
};


struct TLP_RANDGEN{
	
    const TLP_ALGORITHM *ALGINFO;
    void *CURSTATE;
    
};


void TLP_ALGORITHM_SET (TLP_ALGORITHM* ALGINFO, int ALG)
{
	//------------------------------------------------------------
	//	CURRENTLY POSSIBLE ALGORITHMS ARE            
	//
	//	0. 1st GEN. RANLUX223 
	//	1. 1st GEN. RANLUX389
	//	2. 2nd GEN. DOUBLE PRECISION RANLUXD202
	//	3. 2nd GEN. DOUBLE PRECISION RANLUXD397 (HIGHEST LEVEL OF LUXURY)
	//------------------------------------------------------------
	
	if (ALG==0){
		
		ALGINFO->ALGORITHM=0;
		ALGINFO->MAX_INT=0x00ffffffUL;
		ALGINFO->MIN_INT=0;
		ALGINFO->STATE_SIZE=sizeof(RANLUX_STATE);
		ALGINFO->SET_RND=&RANLUX_INIT223;
		ALGINFO->GET_INT=&RANLUX_GEN_INT;
		ALGINFO->GET_DOUBLE=&RANLUX_GEN_DOUBLE;
			
	}else if (ALG==1){
		
		ALGINFO->ALGORITHM=1;
		ALGINFO->MAX_INT=0x00ffffffUL;
		ALGINFO->MIN_INT=0;
		ALGINFO->STATE_SIZE=sizeof(RANLUX_STATE);
		ALGINFO->SET_RND=&RANLUX_INIT389;
		ALGINFO->GET_INT=&RANLUX_GEN_INT;
		ALGINFO->GET_DOUBLE=&RANLUX_GEN_DOUBLE;
		
	}else if (ALG==2){
		
		ALGINFO->ALGORITHM=2;
		ALGINFO->MAX_INT=0xffffffffUL;
		ALGINFO->MIN_INT=0;
		ALGINFO->STATE_SIZE=sizeof(RANLUX_STATE);
		ALGINFO->SET_RND=&RANLUXD_INIT202;
		ALGINFO->GET_INT=&RANLUXD_GEN_INT;
		ALGINFO->GET_DOUBLE=&RANLUXD_GEN_DOUBLE;
		
	}else if (ALG==3){
		
		ALGINFO->ALGORITHM=3;
		ALGINFO->MAX_INT=0xffffffffUL;
		ALGINFO->MIN_INT=0;
		ALGINFO->STATE_SIZE=sizeof(RANLUX_STATE);
		ALGINFO->SET_RND=&RANLUXD_INIT397;
		ALGINFO->GET_INT=&RANLUXD_GEN_INT;
		ALGINFO->GET_DOUBLE=&RANLUXD_GEN_DOUBLE;
		
	}else{
		
		fprintf(stderr, "RANDOM NUMBER GENERATOR ALGORITHM %d NOT IMPLEMENTED\n", ALG); exit(-1);
		
	}
	
}

void TLP_RAND_SET_SEED (TLP_RANDGEN * RANDGEN, unsigned long int EXTSEED)
{
	
   RANDGEN->ALGINFO->SET_RND(RANDGEN->CURSTATE, EXTSEED);
}

TLP_RANDGEN* TLP_RANDGEN_INITIALIZE( TLP_ALGORITHM *ALGINFO, unsigned long int SEED )
{
	
	TLP_RANDGEN *TMP = (TLP_RANDGEN*) malloc (sizeof(TLP_RANDGEN));
		
	TMP->CURSTATE = malloc( ALGINFO->STATE_SIZE );

	TMP->ALGINFO= ALGINFO;
			
	TLP_RAND_SET_SEED ( TMP, SEED);
	
	return TMP;
	
}


inline unsigned long int TLP_RAND_GEN_INT (const TLP_RANDGEN * RANDGEN);
inline double TLP_RAND_GEN_DOUBLE_UNIFORM(const TLP_RANDGEN * RANDGEN);
inline unsigned long int TLP_RAND_GEN_INT_UNIFORM (const TLP_RANDGEN * RANDGEN, unsigned long int MAXINT);

inline unsigned long int TLP_RAND_GEN_INT (const TLP_RANDGEN * RANDGEN)
{
  return (RANDGEN->ALGINFO->GET_INT) (RANDGEN->CURSTATE);
}

inline double TLP_RAND_GEN_DOUBLE_UNIFORM (const TLP_RANDGEN * RANDGEN)
{
  return (RANDGEN->ALGINFO->GET_DOUBLE) (RANDGEN->CURSTATE);
}

inline unsigned long int TLP_RAND_GEN_INT_UNIFORM (const TLP_RANDGEN * RANDGEN, unsigned long int MAXINT)
{
  unsigned long int offset = RANDGEN-> ALGINFO-> MIN_INT;
  unsigned long int range = RANDGEN->ALGINFO->MAX_INT - offset;
  unsigned long int scale;
  unsigned long int k;

  if (MAXINT > range || MAXINT == 0) 
    {
     fprintf(stderr, "CHOICE OF MAXINTEGER %ld LARGER THAN RANGE OF GENERATOR %ld or 0\n", MAXINT,range); exit(-1);
    }

  scale = range / MAXINT;

  do
    {
      k = (((RANDGEN->ALGINFO->GET_INT) (RANDGEN->CURSTATE)) - offset) / scale;
    }
  while (k >= MAXINT);

  return k;
}



//	************************************************************************************************************************** //
//	****                                                                                                                  **** //
//	****                        THE GAUSS DISTRIBUTION ACCORDING TO GSL IMPLEMENTATION OF GAUSS ZIGGURAT                  **** //
//	****                                                                                                                  **** //
//	************************************************************************************************************************** //

/* position of right-most step */
#define PARAM_R 3.44428647676

/* tabulated values for the heigt of the Ziggurat levels */
static const double ytab[128] = {
  1, 0.963598623011, 0.936280813353, 0.913041104253,
  0.892278506696, 0.873239356919, 0.855496407634, 0.838778928349,
  0.822902083699, 0.807732738234, 0.793171045519, 0.779139726505,
  0.765577436082, 0.752434456248, 0.739669787677, 0.727249120285,
  0.715143377413, 0.703327646455, 0.691780377035, 0.68048276891,
  0.669418297233, 0.65857233912, 0.647931876189, 0.637485254896,
  0.62722199145, 0.617132611532, 0.607208517467, 0.597441877296,
  0.587825531465, 0.578352913803, 0.569017984198, 0.559815170911,
  0.550739320877, 0.541785656682, 0.532949739145, 0.524227434628,
  0.515614886373, 0.507108489253, 0.498704867478, 0.490400854812,
  0.482193476986, 0.47407993601, 0.466057596125, 0.458123971214,
  0.450276713467, 0.442513603171, 0.434832539473, 0.427231532022,
  0.419708693379, 0.41226223212, 0.404890446548, 0.397591718955,
  0.390364510382, 0.383207355816, 0.376118859788, 0.369097692334,
  0.362142585282, 0.355252328834, 0.348425768415, 0.341661801776,
  0.334959376311, 0.328317486588, 0.321735172063, 0.31521151497,
  0.308745638367, 0.302336704338, 0.29598391232, 0.289686497571,
  0.283443729739, 0.27725491156, 0.271119377649, 0.265036493387,
  0.259005653912, 0.253026283183, 0.247097833139, 0.241219782932,
  0.235391638239, 0.229612930649, 0.223883217122, 0.218202079518,
  0.212569124201, 0.206983981709, 0.201446306496, 0.195955776745,
  0.190512094256, 0.185114984406, 0.179764196185, 0.174459502324,
  0.169200699492, 0.1639876086, 0.158820075195, 0.153697969964,
  0.148621189348, 0.143589656295, 0.138603321143, 0.133662162669,
  0.128766189309, 0.123915440582, 0.119109988745, 0.114349940703,
  0.10963544023, 0.104966670533, 0.100343857232, 0.0957672718266,
  0.0912372357329, 0.0867541250127, 0.082318375932, 0.0779304915295,
  0.0735910494266, 0.0693007111742, 0.065060233529, 0.0608704821745,
  0.056732448584, 0.05264727098, 0.0486162607163, 0.0446409359769,
  0.0407230655415, 0.0368647267386, 0.0330683839378, 0.0293369977411,
  0.0256741818288, 0.0220844372634, 0.0185735200577, 0.0151490552854,
  0.0118216532614, 0.00860719483079, 0.00553245272614, 0.00265435214565
};

/* tabulated values for 2^24 times x[i]/x[i+1],
 * used to accept for U*x[i+1]<=x[i] without any floating point operations */
static const unsigned long ktab[128] = {
  0, 12590644, 14272653, 14988939,
  15384584, 15635009, 15807561, 15933577,
  16029594, 16105155, 16166147, 16216399,
  16258508, 16294295, 16325078, 16351831,
  16375291, 16396026, 16414479, 16431002,
  16445880, 16459343, 16471578, 16482744,
  16492970, 16502368, 16511031, 16519039,
  16526459, 16533352, 16539769, 16545755,
  16551348, 16556584, 16561493, 16566101,
  16570433, 16574511, 16578353, 16581977,
  16585398, 16588629, 16591685, 16594575,
  16597311, 16599901, 16602354, 16604679,
  16606881, 16608968, 16610945, 16612818,
  16614592, 16616272, 16617861, 16619363,
  16620782, 16622121, 16623383, 16624570,
  16625685, 16626730, 16627708, 16628619,
  16629465, 16630248, 16630969, 16631628,
  16632228, 16632768, 16633248, 16633671,
  16634034, 16634340, 16634586, 16634774,
  16634903, 16634972, 16634980, 16634926,
  16634810, 16634628, 16634381, 16634066,
  16633680, 16633222, 16632688, 16632075,
  16631380, 16630598, 16629726, 16628757,
  16627686, 16626507, 16625212, 16623794,
  16622243, 16620548, 16618698, 16616679,
  16614476, 16612071, 16609444, 16606571,
  16603425, 16599973, 16596178, 16591995,
  16587369, 16582237, 16576520, 16570120,
  16562917, 16554758, 16545450, 16534739,
  16522287, 16507638, 16490152, 16468907,
  16442518, 16408804, 16364095, 16301683,
  16207738, 16047994, 15704248, 15472926
};

/* tabulated values of 2^{-24}*x[i] */
static const double wtab[128] = {
  1.62318314817e-08, 2.16291505214e-08, 2.54246305087e-08, 2.84579525938e-08,
  3.10340022482e-08, 3.33011726243e-08, 3.53439060345e-08, 3.72152672658e-08,
  3.8950989572e-08, 4.05763964764e-08, 4.21101548915e-08, 4.35664624904e-08,
  4.49563968336e-08, 4.62887864029e-08, 4.75707945735e-08, 4.88083237257e-08,
  5.00063025384e-08, 5.11688950428e-08, 5.22996558616e-08, 5.34016475624e-08,
  5.44775307871e-08, 5.55296344581e-08, 5.65600111659e-08, 5.75704813695e-08,
  5.85626690412e-08, 5.95380306862e-08, 6.04978791776e-08, 6.14434034901e-08,
  6.23756851626e-08, 6.32957121259e-08, 6.42043903937e-08, 6.51025540077e-08,
  6.59909735447e-08, 6.68703634341e-08, 6.77413882848e-08, 6.8604668381e-08,
  6.94607844804e-08, 7.03102820203e-08, 7.11536748229e-08, 7.1991448372e-08,
  7.2824062723e-08, 7.36519550992e-08, 7.44755422158e-08, 7.52952223703e-08,
  7.61113773308e-08, 7.69243740467e-08, 7.77345662086e-08, 7.85422956743e-08,
  7.93478937793e-08, 8.01516825471e-08, 8.09539758128e-08, 8.17550802699e-08,
  8.25552964535e-08, 8.33549196661e-08, 8.41542408569e-08, 8.49535474601e-08,
  8.57531242006e-08, 8.65532538723e-08, 8.73542180955e-08, 8.8156298059e-08,
  8.89597752521e-08, 8.97649321908e-08, 9.05720531451e-08, 9.138142487e-08,
  9.21933373471e-08, 9.30080845407e-08, 9.38259651738e-08, 9.46472835298e-08,
  9.54723502847e-08, 9.63014833769e-08, 9.71350089201e-08, 9.79732621669e-08,
  9.88165885297e-08, 9.96653446693e-08, 1.00519899658e-07, 1.0138063623e-07,
  1.02247952126e-07, 1.03122261554e-07, 1.04003996769e-07, 1.04893609795e-07,
  1.05791574313e-07, 1.06698387725e-07, 1.07614573423e-07, 1.08540683296e-07,
  1.09477300508e-07, 1.1042504257e-07, 1.11384564771e-07, 1.12356564007e-07,
  1.13341783071e-07, 1.14341015475e-07, 1.15355110887e-07, 1.16384981291e-07,
  1.17431607977e-07, 1.18496049514e-07, 1.19579450872e-07, 1.20683053909e-07,
  1.21808209468e-07, 1.2295639141e-07, 1.24129212952e-07, 1.25328445797e-07,
  1.26556042658e-07, 1.27814163916e-07, 1.29105209375e-07, 1.30431856341e-07,
  1.31797105598e-07, 1.3320433736e-07, 1.34657379914e-07, 1.36160594606e-07,
  1.37718982103e-07, 1.39338316679e-07, 1.41025317971e-07, 1.42787873535e-07,
  1.44635331499e-07, 1.4657889173e-07, 1.48632138436e-07, 1.50811780719e-07,
  1.53138707402e-07, 1.55639532047e-07, 1.58348931426e-07, 1.61313325908e-07,
  1.64596952856e-07, 1.68292495203e-07, 1.72541128694e-07, 1.77574279496e-07,
  1.83813550477e-07, 1.92166040885e-07, 2.05295471952e-07, 2.22600839893e-07
};


double TLP_RAND_GAUSSIAN_ZIGGURAT(const TLP_RANDGEN * RANDGEN, const double sigma)
{
  unsigned long int i, j;
  int sign;
  double x, y;

  const unsigned long int range = RANDGEN->ALGINFO->MAX_INT - RANDGEN->ALGINFO->MIN_INT;
  const unsigned long int offset = RANDGEN->ALGINFO->MIN_INT;

  while (1)
    {
      if (range >= 0xFFFFFFFF)
        {
          unsigned long int k = TLP_RAND_GEN_INT(RANDGEN) - offset;
          i = (k & 0xFF);
          j = (k >> 8) & 0xFFFFFF;
        }
      else if (range >= 0x00FFFFFF)
        {
          unsigned long int k1 = TLP_RAND_GEN_INT(RANDGEN) - offset;
          unsigned long int k2 = TLP_RAND_GEN_INT(RANDGEN) - offset;
          i = (k1 & 0xFF);
          j = (k2 & 0x00FFFFFF);
        }
      else
        {
          i = TLP_RAND_GEN_INT_UNIFORM (RANDGEN, 256); /*  choose the step */
          j = TLP_RAND_GEN_INT_UNIFORM (RANDGEN, 16777216);  /* sample from 2^24 */
        }

      sign = (i & 0x80) ? +1 : -1;
      i &= 0x7f;

      x = j * wtab[i];

      if (j < ktab[i])
        break;

      if (i < 127)
        {
          double y0, y1, U1;
          y0 = ytab[i];
          y1 = ytab[i + 1];
          U1 = TLP_RAND_GEN_DOUBLE_UNIFORM (RANDGEN);
          y = y1 + (y0 - y1) * U1;
        }
      else
        {
          double U1, U2;
          U1 = 1.0 - TLP_RAND_GEN_DOUBLE_UNIFORM (RANDGEN);
          U2 = TLP_RAND_GEN_DOUBLE_UNIFORM (RANDGEN);
          x = PARAM_R - log (U1) / PARAM_R;
          y = exp (-PARAM_R * (x - 0.5 * PARAM_R)) * U2;
        }

      if (y < exp (-0.5 * x * x))
        break;
    }

  return sign * sigma * x;
}


