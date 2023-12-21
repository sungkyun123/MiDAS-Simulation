
using namespace std;

struct REQUEST{
    double time;
    int type;
    int lba;
    int io_size;
    int stream;
    double timestamp;
    bool sw;
};

double age_function(struct SSD *ssd, struct STATS *stats, int seg_idx);

int CB_VS(struct SSD *ssd, struct STATS *stats, class GROUP *group[], int gc_group);
int FIFO_VS(class GROUP *group[], int gc_group);
int GREEDY_VS(struct SSD *ssd, class GROUP *group[], int gc_group);
int victim_selection(struct SSD *ssd, struct STATS *stats, class GROUP *group[], int gc_group);
void update_gc_results(struct SSD* ssd, struct STATS* stats, int local_copy, int victim_gnum, int victim_seg);
void initialize_segment(struct SSD *ssd, int victim_idx, int victim_seg, int victim_gnum);
int do_gc(struct SSD *ssd, struct STATS *stats, class GROUP *group[], int gc_group);
int gc_victim_group(struct SSD *ssd, class GROUP *group[]);
int GC(struct SSD* ssd, struct STATS *stats, class GROUP *group[]);

int write(int lba, struct SSD *ssd, struct STATS *stats, class GROUP *group[], bool sw);
int trim(int lba, struct SSD *ssd, struct STATS *stats);


void req_processing(char* raw_req, struct REQUEST *req, struct STATS *stats);
int submit_io(struct SSD *ssd, struct STATS *stats, class GROUP* group[], struct REQUEST *req);

int stats_clear(struct SSD* ssd, struct STATS* stats);
void stat_update(struct STATS* stats, int ret, int type);
int simul_info(struct SSD *ssd, struct STATS *stats, class GROUP* group[], int time_gap);
int display_result(struct SSD *ssd, struct STATS *stats, char *workload, char *vs_policy);
int display_debug(struct SSD *ssd, char *workload);
int logging_result(struct SSD *ssd, struct STATS *stats, char *workload, char *vs_policy);


int ssd_simulation_sw(struct SSD* ssd, struct STATS* stats, class GROUP* group[], char *workload, int group_size[]);
int ssd_simulation(struct SSD* ssd, struct STATS* stats, class GROUP* group[], char* workload, int group_size[]);




