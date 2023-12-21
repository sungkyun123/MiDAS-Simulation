#include <pthread.h>
#include <iostream>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <queue>
#include <deque>

using namespace std;

#define WRITE 0
#define REMOVE 1

typedef struct miniature_model mini_model;

/* manage time window and interval unit size
 * time: # of interval count
 */
typedef struct model_time {
        bool is_real;
        //unsigned long long time_window;
        unsigned long long current_time; //time, # of interval count
        uint32_t request_time; //time, # of requests, max: 128*512
        unsigned long long load_timing; //to remove sequential write in model
        uint32_t interval_unit_size;
        unsigned long long extra_load; //in real-world workload, wait extra time after LOAD_END signal
}mtime;


/* used when modeling is over & get new group configuration
 * valid: true if the configuration information is valid
 * when midas use that info and change confi, valid must turn off to false
 */
typedef struct group_configuration{
        bool valid; // true if valid info
	int model_ver; //TODO remove
	int gnum;
        int M_gnum;
	double *vr;
        double *commit_vr; 
        int commit_g;
        uint32_t* gsize;
        double WAF;
	bool err_stat;
	double g0_traffic;
}G_INFO;


struct miniature_model {
	uint32_t lba_sampling_ratio; //sampling ratio
	uint32_t interval_unit_size; //1 interval unit size, # of segment
	unsigned long long desired_tw;
	unsigned long long time_window; // time window size, # of interval count
	unsigned long long real_tw; //time window size, # of pages
	uint32_t entry_num;

	uint32_t fnumber; //number of checking first interval
	uint32_t* checking_first_interval; //interval counts to check first interval
	uint32_t first_count; //update counts for checked first interval
	
	unsigned long long first_interval; //for new first interval analyzer
	pthread_t thread_id;

	unsigned long long* time_stamp; //time stamp
	unsigned long long* model_count; // counting update count per intervals
	uint32_t* updated_lbas; // for first interval's LBAs

	int one_op;
	unsigned long long hot_req_count;
	unsigned long long tot_req_count;
	double hot_lba_count;
	double tot_lba_count;
	int modeling_num=-1;

	long no_access_lba;
	unsigned long long total_count;
	bool model_on;
	mtime *mm_time;
	int model_idx;
	double (*waf_predictor) (mini_model *mmodel, double *valid_ratio_list, int group);

	int err_cnt;
	int tmp_err_cnt;
};

class MODEL_Q{
	public:
		deque<double> g0_traffic_queue;
		deque<double> g0_size_queue;
		deque<double> g0_valid_queue;

		double g0_traffic;
		double g0_size;
		double g0_valid;
		int queue_max;

		int extra_size;
		double extra_traffic;

		double calc_traffic;
		int calc_size;
		int calc_unit;

		int best_extra_size;
		double best_extra_traffic;
		double best_extra_unit;
};

double dqueue_pop_back(deque<double> *dqueue);
void model_create(mini_model *mmodel, int write_size, int type);
void time_managing(mini_model *mmodel, char mode);
void model_initialize(mini_model *mmodel, bool on);
void model_resizing_timewindow(mini_model *mmodel, unsigned long long timewindow);
int check_time_window(mini_model *mmodel, char mode);
int check_interval(mini_model *mmodel, uint32_t lba, char mode);
int check_first_interval(mini_model *mmodel, uint32_t lba, char mode);
void *making_group_configuration(void *arg);
void *making_group_configuration2(void *arg);
void *making_group_configuration3(void *arg);
void print_config(int, uint32_t*, double, double*);
void print_config_new(int, uint32_t*, double, double*, double);
void print_config_into_log(int, uint32_t*, double, double*);
void print_config2(int, uint32_t*, double, double*);
void print_config3(int, uint32_t*, double, double*, double);
double WAF_predictor(mini_model *mmodel, double *, int);
double WAF_predictor_2(mini_model *mmodel, double *, int);
double one_WAF_predictor(mini_model *mmodel, double *, int);
double two_WAF_predictor(mini_model *mmodel, double *, int);
unsigned long long resizing_model(mini_model *mmodel);
double *one_valid_ratio_predictor(mini_model *mmodel, uint32_t*, uint32_t, unsigned long long, int);
double *two_valid_ratio_predictor(mini_model *mmodel, uint32_t*, uint32_t, unsigned long long, int);
double *valid_ratio_predictor(mini_model *mmodel, uint32_t*, uint32_t, unsigned long long);
void initialize_first_interval(mini_model *mmodel);
void remove_first_interval(mini_model *mmodel);
void model_destroy(mini_model *mmodel);
int update_count(mini_model *mmodel, uint32_t lba, char mode);

// Desnoyer

double FIFO_predict(double op);
double None_FIFO_predict(mini_model *mmodel, double op);
double one_group_predictor(mini_model *mmodel, uint32_t*, uint32_t, unsigned long long);
