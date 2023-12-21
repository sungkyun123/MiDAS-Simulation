#include <string>
#include <sstream>
#include <time.h>
#include <cmath>
#include <signal.h>
#include <deque>
#include "model.h"
#include "ssd_config.h"
#include "queue.h"
#include <cstring>
#include "hf.h"

extern SSD_SPEC *ssd_spec;
extern struct SSD *ssd;
extern mini_model *mmodel;
extern MODEL_Q *model_q;
extern struct STATS *stats;

//Hot filter data structure inititialization
void hf_init(struct SSD* ssd, struct HotFilter **hotf){
        (*hotf) = (struct HotFilter*)malloc(sizeof(struct HotFilter));
        (*hotf)->make_flag = 0;
        (*hotf)->use_flag = 0;
        (*hotf)->ready_flag=0;

        //TODO check max value of hotfilter
        //(*hotf)->max_val = 1; // 2-bit per LBA
        (*hotf)->max_val = 3;
        (*hotf)->hot_val=3;
        (*hotf)->tw_ratio = 1.0;

        (*hotf)->tw = (long)ssd_spec->LOGSIZE/(long)4096;

        (*hotf)->left_tw = (*hotf)->tw;
        (*hotf)->cold_tw = (long)ssd_spec->LOGSIZE/(long)4096;
        //(*hotf)->cold_tw = 0;

        (*hotf)->G0_vr_num = 0;
        (*hotf)->G0_vr_sum = 0;
        (*hotf)->seg_age=0.0;
        (*hotf)->seg_num=0.0;
        (*hotf)->avg_seg_age=0.0;

        (*hotf)->G0_traffic_ratio = 0;
        (*hotf)->G0_traffic = 0;
        (*hotf)->tot_traffic = 0;
        (*hotf)->hot_lba_num = 0;
        (*hotf)->valid_lba_num=0.;

        (*hotf)->cur_hf = (int*)malloc(sizeof(int)*ssd_spec->LBANUM);
        //hotf->new_hf = (int*)malloc(sizeof(int)*ssd_spec->LBANUM);

        memset((*hotf)->cur_hf, 0, sizeof(int)*ssd_spec->LBANUM);
        //memset(hotf->new_hf, 0, sizeof(int)*ssd_spec->LBANUM);

        (*hotf)->err_cnt=0;
}

void hf_destroy(struct HotFilter *hotf) {
        free(hotf->cur_hf);
        //free(hotf->new_nf);
        free(hotf);
}

//Metadata Reset for next hot filter generation
void hf_metadata_reset(struct HotFilter *hotf){
        hotf->use_flag = 0;
        hotf->make_flag = 0;
        hotf->ready_flag=0;
        hotf->G0_vr_num = 0;
        hotf->G0_vr_sum = 0;
        hotf->seg_age=0.0;
        hotf->seg_num=0.0;
        hotf->avg_seg_age=0.0;
        hotf->G0_traffic_ratio = 0;
        hotf->G0_traffic = 0;
        hotf->tot_traffic = 0;
        hotf->hot_lba_num=0;
        hotf->left_tw = hotf->tw;
}

//Hot filter reset default: flag->0
void hf_reset(int flag, struct HotFilter *hotf){
        memset(hotf->cur_hf, 0, sizeof(int)*ssd_spec->LBANUM);
        hotf->left_tw = hotf->tw;        
        if (flag == 1) hotf->make_flag = 1; // not used 
        hotf->use_flag = 0;
}

void hf_update_model(double traffic, struct HotFilter *hotf) {
        if (mmodel->model_on == false) return;
        double tmp;
        if (traffic == 0.0) return;
        if ((int)model_q->g0_traffic_queue.size() >= model_q->queue_max) {
                tmp = dqueue_pop_back(&model_q->g0_traffic_queue);
        }
        if ((int)model_q->g0_valid_queue.size() >= model_q->queue_max) {
                tmp = dqueue_pop_back(&model_q->g0_valid_queue);
        }
        if ((int)model_q->g0_size_queue.size() >= model_q->queue_max) {
                tmp = dqueue_pop_back(&model_q->g0_size_queue);
        }

        model_q->g0_traffic_queue.push_front(traffic);
        model_q->g0_valid_queue.push_front(hotf->G0_vr_sum/hotf->G0_vr_num);
        model_q->g0_size_queue.push_front(ssd->TOTAL_GNUM[0]);

        return;
}

void hf_convert_generate_to_run(struct HotFilter *hotf) {
        //change flags

        if (hotf->left_tw < 0) {
                hotf->left_tw=0;
        }

        //print HotFilter info
        double prev_seg_age=0.0;
        double tr=0.0;
        tr = hotf->G0_traffic/hotf->tot_traffic;
        hotf->G0_traffic_ratio=tr;

        if (hotf->G0_traffic != 0) hf_update_model(tr, hotf);
        prev_seg_age = hotf->tw/hotf->tw_ratio;

        //g0_size: G0 size ratio of total segment size
        //avg_g0 : average age of G0 victim segment
        double avg_g0 = 0.0;
        if (hotf->seg_num==0) {
                printf("there is no G0 gc segment (in hf generate())\n");
                abort();
        }
        hotf->avg_seg_age = hotf->seg_age/hotf->seg_num;
        avg_g0 = hotf->avg_seg_age/(double)ssd_spec->PPS;

        //resize the time window of new hotfilter
        //G0 average age * tw ratio (it is related to the max bit of hotfilter)
        //TODO check avg_seg_age update
        //2

        if (prev_seg_age != 0.0) {
                if (hotf->avg_seg_age >= prev_seg_age*2) {
                        //limits the new avg_seg_age
                        hotf->avg_seg_age = prev_seg_age*1.5;
                } else {
                        //calculate the average age
                        hotf->avg_seg_age = (hotf->avg_seg_age+prev_seg_age)/2.0;
                }
        }


        //if (prev_seg_age != 0.0) hotf->avg_seg_age = (hotf->avg_seg_age+prev_seg_age)/2.0;
        avg_g0 = hotf->avg_seg_age/(double)ssd_spec->PPS;
        long tw = (long)(hotf->avg_seg_age*hotf->tw_ratio);

        hotf->tw = tw;
        hotf->left_tw = tw-1;
        hotf->seg_age=0.0;
        hotf->seg_num=0.0;
        hotf->G0_traffic=0.0;
        hotf->tot_traffic=0.0;
        hotf->G0_vr_sum=0;
        hotf->G0_vr_num=0;

        return;
}

//Hot filter generate code
//old seg: segment that the lba was placed before
//user_group: group number that the lba is moving
void hf_generate(struct SSD *ssd, int lba, int old_seg, class GROUP *group[], struct HotFilter *hotf, int hflag){
        //if the hotfilter is not generating
        if (hotf->make_flag == 0) {
                return;
        }
        if (hotf->cold_tw > 0) {
                hotf->cold_tw--;
                return;
        }
        //generating hotfilter is over
        //TODO check converting hotfilter info
        //1
        if ((hotf->make_flag == 1) && ((hotf->left_tw<=0) && (hotf->seg_num > 0))){ //TW for hot filter generation is finished
        //if ((hotf->make_flag==1) && (hotf->seg_num > 15)) {
                hf_convert_generate_to_run(hotf);
                //return;
        }
        if (hotf->make_flag == 1 && hotf->left_tw == hotf->tw) { //generation start condition
                //printf("[HF-NOTICE] HOT FILTER CREATION START\n");
        }
        //update the LBA to the hotfilter
        if (hflag) {
                hotf->left_tw --; //decrease counter of tw
                //if the tw is over, but there is no erased segment in group 0
                //TODO decide to increase the tw
                //if (hotf->left_tw < 0) hotf->tw++;
        }


        if (hflag) {
                if (ssd->gnum_info[old_seg] == 0) {
                        //LBA in group 0 is pretended HOT regardless of the age
                        if (hotf->cur_hf[lba] <= hotf->max_val-1) {
                                hotf->cur_hf[lba]++; //increase counter of the LBA
                                if (hotf->cur_hf[lba] == hotf->hot_val) hotf->hot_lba_num++;
                        }
                } else if (ssd->gnum_info[old_seg] == 1) {
                        //LBA in group 1 is not in current hotfilter
                        //but can be add to new hotfilter, if the age is lower than the G0 age
                        double seg_age = 0.0;
                        if (ssd->seg_stamp[old_seg] == -1.0) {
                                //open segment
                                seg_age = 0.0;
                        }
                        else seg_age = (double)(stats->cur_wp - ssd->seg_stamp[old_seg]);
                        //(hotf->tw / tw_ratio) is the average G0 age
                        double tmp_age = hotf->tw/hotf->tw_ratio;
                        if (seg_age <= tmp_age) {
                                if (hotf->cur_hf[lba] <= hotf->max_val-1) {
                                        hotf->cur_hf[lba]++;
                                        if (hotf->cur_hf[lba] == hotf->hot_val) hotf->hot_lba_num++;
                                }
                        }
                        else {
                                //invalid lba in group 1 victim, but the age is too long
                                if (hotf->cur_hf[lba] > 0) {
                                        hotf->cur_hf[lba]--;
                                        //if (hotf->cur_hf[lba] == (hotf->max_val-1)) hotf->hot_lba_num--;
                                        if (hotf->cur_hf[lba] == (hotf->hot_val-1)) hotf->hot_lba_num--;
                                }
                        }
                }
        } else {
                //valid lba in group 0 victim
                //TODO --? or avg_age check?
                if (hotf->cur_hf[lba] > 0) {
                        hotf->cur_hf[lba]--;
                        if (hotf->cur_hf[lba] == hotf->hot_val-1) hotf->hot_lba_num--;
                        if (ssd->gnum_info[old_seg] == 0) hotf->cur_hf[lba]--; //hot penalty
                }
        }
}

void hf_calculate_valid_lba(int lba, int group_num, char vflag, struct HotFilter *hotf) {
        if (vflag == 0) {
                //write (valid)
                if (group_num == 0) hotf->valid_lba_num++;
        } else if (vflag==1) {
                //write (invalid)
                if (group_num == 0) hotf->valid_lba_num--;
        } else {
                //gc (removed)
                hotf->valid_lba_num--;
        }
        return;
}

//check wthether the LBA is assigned to G0 or G1
int hf_check(int lba, struct HotFilter *hotf){
        if (hotf->use_flag == 0) return 1; //if hot filter not available
        else { //if hot filter available
                hotf->tot_traffic ++; //increase total traffic
                if (hotf->cur_hf[lba] > hotf->max_val) {
                        printf("hotfilter value over!!! (LBA: %d, value: %d)\n", lba, hotf->cur_hf[lba]);
                        abort();
                } else if (hotf->cur_hf[lba] < 0) {
                        printf("hotfilter value over!!! (LBA: %d, value: %d)\n", lba, hotf->cur_hf[lba]);
                        abort();
                }
                if (hotf->cur_hf[lba] >= hotf->hot_val) { //if the LBA is hot
                        hotf->G0_traffic ++; //increase the traffic of G0
                        return 0; //HOT
                }
                else return 1; //G1
        }
}


