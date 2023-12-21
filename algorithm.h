using namespace std;

int do_gc_for_split(struct SSD *ssd, struct STATS *stats, class GROUP *group[], int gc_group);
int error_cmp(struct SSD *ssd, class GROUP *group[]);
int group_split(struct SSD *ssd, struct STATS *stats, class GROUP *group[], int new_last);
void active_split(struct SSD *ssd, struct STATS *stats, class GROUP *group[]);
int mida_off(struct SSD *ssd, struct STATS *stats, class GROUP *group[]);
int active_merge(struct SSD *ssd, struct STATS *stats, class GROUP *group[], int target_g, int victim_g);
void hot_merge(struct SSD* ssd, struct STATS* stats, class GROUP *group[]);
void merge_group(struct SSD *ssd, struct STATS *stats, class GROUP *group[]);
int gsize_check(struct SSD* ssd, class GROUP *group[]);
int move_margin_seg(struct SSD* ssd, struct STATS *stats, class GROUP *group[], int victim_g, int target_g, int size);
int change_config_2(struct SSD* ssd, struct STATS *stats, class GROUP *group[]);
int checking_ginfo(struct SSD* ssd, struct STATS* stats, class GROUP *group[]);
void checking_err_time(struct SSD* ssd, struct STATS* stats, class GROUP *group[]);
void modeling_check(struct SSD* ssd, struct STATS* stats, class GROUP *group[], int user_group, int lba);
void checking_config_apply(struct SSD* ssd, struct STATS* stats, class GROUP *group[], int gc_group);
int infinite_gc_handling(struct SSD* ssd, struct STATS* stats, class GROUP *group[], int gc_group);
void group_number_check(struct SSD *ssd);
int check_applying_config(struct SSD *ssd, struct STATS *stats, class GROUP* group[]);
int mida_on_off(struct SSD *ssd, struct STATS *stats, class GROUP* group[], int SWITCH);
void err_handling(struct SSD* ssd, struct STATS* stats, class GROUP* group[]);

