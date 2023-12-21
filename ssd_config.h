using namespace std;

struct SSD {
	bool *itable; // valid: 0, invalid: 1
	int *mtable; // mapping table (idx: lba, value: physical position)
	int *irtable; // invalid page count of each block
	unsigned char *gnum_info; // group number of each block
	bool *fill_info; // filled: 0, not filled: 1
	double *VPG; // valid ratio per group
	double *EPG; // erase count per group
	int **active; // active block per group
	double *TMP_VPG; // valid ratio per group with time window
	double *TMP_EPG; // erase count with time window
	double *DEBUG_VPG;
	double *DEBUG_EPG;
	double *CUR_VPG; // valid ratio per group with some state
	double *TOTAL_GNUM; // # of blocks in each group
	double *seg_stamp; // time per block (not used yet)
	struct OOB *oob; // oob per page
	double *err_VPG;
	double *err_EPG;
	long **lifespan;
	long *lba_life;
	int config_version; //0: MiDASS, 1: MiDASS+HF
	int gnum;
	
};

struct SSD_SPEC {
	int vs_policy; //victim slection policy 0: default (MiDAS)
	int GROUPNUM; //# of groups
	unsigned long long int LOGSIZE; //logical size of SSD
	unsigned long long int DEVSIZE; //physical size of SSD
	int dev_gb;
	int seg_mb;
	double OP; //overprovisioning rate
	int PPB; //pages per block
	int BPS; //blocks per segment
	int PPS; //pages per segment
	int PGSIZE; // page size (4KB)
	int BLKNUM; // # of blocks in SSD
	int BLKSIZE; // block size (PPB*PGSIZE)
	int PGNUM; // BLKNUM*PPB
	int SEGSIZE; // segment size (BPS*BLKSIZE)
	int SEGNUM; // # of segments in SSD
	int fifo_mode; // victim selection of last group
	int LBANUM; // # of LBAs in SSD (LOGSIZE/PGSIZE)
	int FREENUM; // # of free segments in SSD
	int MAXGNUM; // Maximum # of groups for MiDAS
	int naive_start;
	int mida_on;
	int INF_GC_CNT;
	int NAIVE_N;
	
};

struct STATS {
	unsigned long long cur_req; //current # of request
	double cur_time; //current time in SSD
	unsigned long long cur_wp; //current # of write page
	unsigned long long load_wp;
	unsigned long long tot_req; //total # of request
	unsigned long long copy; //current # of copy page
	int tmp_wp; //current # of write page in time window
	int tmp_copy; //current # of copy page in time window
	int erase; //current # of block erase
	int progress;
	unsigned long long vp;
	double util;
	double fourty_waf;
	double fourty_wp;
	long m_wp;
	long m_copy;
	unsigned long long unknown_req;
	unsigned long long unknown_cnt;
	double checktime;
	int tw;
	long err_write;
	int err_flag;

	int commit_bit;
	double tmp_waf;
};

struct OOB {
	int lba; //lba of page
	int gnum; //group number of page
};

class GROUP{
	public:
		deque<int> fill_queue; //group management queue
		int size; //group size
};

void ssd_init(struct SSD* ssd, class GROUP* group[], int gnum, int vs_policy, int naive_start, int dev_gb, int seg_mb);
void group_init(struct SSD* ssd, class GROUP *group[], int gnum, int mida_on, int size[]);
void stats_init(struct STATS *stats);
int policy_to_int(char *policy);

