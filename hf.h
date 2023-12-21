using namespace std;

typedef struct HotFilter {
        int *cur_hf; // hot filter
        int max_val; //maximum value of bits per LBA (if 1bit: 1, 2bit: 3, 3bit: 7)
        int hot_val;
        double tw_ratio;
        int make_flag; //generate flag, 0: don't modifiying, 1: modifiying...
        int use_flag; //usage flag, 0: don't used, 1: use!
        int ready_flag;

        long tw; //unit: page
        long left_tw; // left time window;
        long cold_tw;

        double G0_vr_sum;
        double G0_vr_num;

        double seg_age;
        double seg_num;
        double avg_seg_age;

        double G0_traffic_ratio; //G0 traffic ratio for MC model
        double tot_traffic; //total traffic in TW
        double G0_traffic; //G0 traffic in TW
        int hot_lba_num; //number of hot lba (count > max_val)
        double valid_lba_num; //number of valid lba in group 0

        int err_cnt;
 	int hf_cnt;
 	int tmp_err_cnt;
}HF;

void hf_reset(int flag, struct HotFilter *hotf);
void hf_generate(struct SSD *ssd, int lba, int old_seg, class GROUP *group[], struct HotFilter *hotf, int hflag);
void hf_calculate_valid_lba(int lba, int group_num, char vflag, struct HotFilter *hotf);

void hf_init(struct SSD* ssd, struct HotFilter **hotf);
void hf_destroy(struct HotFilter *hotf);
void hf_metadata_reset(struct HotFilter *hotf);

void hf_convert_generate_to_run(struct HotFilter *hotf);
void hf_update_model(double traffic, struct HotFilter *hotf);

void hf_calculate_valid_lba(int lba, int group_num, char vflag, struct HotFilter *hotf);
int hf_check(int lba, struct HotFilter *hotf);

