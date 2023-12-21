#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <stdbool.h>
#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <cmath>
#include <signal.h>
#include <deque>
#include "model.h"
#include "ssd_config.h"
#include "queue.h"
#include "hf.h"
#include "algorithm.h"
#include "ssdsimul.h"

using namespace std;

#define OFF 0 // naive_mida OFF
#define ON 1 // naive_mida ON

#define ILOOP 25 //When the # of recursive GC is over ILOOP, merge the last two groups

extern struct SSD_SPEC *ssd_spec;
extern struct SSD *ssd;
extern struct STATS *stats;

extern deque<int> fbqueue; //global free block management queue

extern mini_model *mmodel;
extern G_INFO *ginfo;
extern MODEL_Q *model_q;

extern struct HotFilter *hf_gen;


void merge_group(struct SSD *ssd, struct STATS *stats, class GROUP *group[]);
int mida_on_off(struct SSD *ssd, struct STATS *stats, class GROUP* group[], int SWITCH);
int gsize_check(struct SSD* ssd, class GROUP *group[]);
int check_applying_config(struct SSD *ssd, struct STATS *stats, class GROUP* group[]);
void err_handling(struct SSD* ssd, struct STATS* stats, class GROUP* group[]);

int do_gc_for_split(struct SSD *ssd, struct STATS *stats, class GROUP *group[], int gc_group){
    int local_copy = 0; //# of copy page in this GC
    int select_group = gc_group;
    if (gc_group >= ssd_spec->naive_start) select_group = ssd_spec->naive_start;
    //whole group is naive MiDA mode (HOT group is naive MiDA too)
    if (ssd_spec->naive_start == 1) select_group = ssd_spec->naive_start;
    int victim_seg = victim_selection(ssd, stats,group, select_group); //victim segment selection
    int victim_idx = victim_seg*ssd_spec->PPS; //victim segment entry point (page)
    int target_group; //destination group of valid page (copy target group)
    int copy_lba; //valid page lba 
    int old_ppa; //ppa of old (invalid) data
    int new_ppa; //new assigned ppa for valid copy
    int victim_gnum = ssd->gnum_info[victim_seg];

    int off_flag = 0;
    if (victim_gnum < ssd_spec->GROUPNUM-1) target_group = victim_gnum + 1; //goto next group if victim group is not last group
    else target_group = victim_gnum; // go back to GC group because the victim group is last group
    //if there is no active segment in target group, then make new active segment
    if (ssd->active[target_group][0] == -1) {
        if (ssd->TOTAL_GNUM[target_group] != 0) {
            //ssd->TOTAL_GNUM[target_group]=0;
            simul_info(ssd, stats, group, 100);
            abort();
        }
        ssd->active[target_group][0] = queue_pop_back(&fbqueue);
        ssd->gnum_info[ssd->active[target_group][0]] = target_group;
        ssd->active[target_group][1] = 0;
    }
    int now_active = ssd->active[target_group][0]; //active segment number of target group


    for (int i = 0; i < ssd_spec->PPS; i++){
        if(ssd->active[target_group][1] == ssd_spec->PPS){ // if active segment is full
            if (ssd->active[target_group][0] == -1) {
                abort();
            }
            if (target_group >= ssd_spec->naive_start) {
                //naive MiDA group
                group[ssd_spec->naive_start]->fill_queue.push_front(now_active);
            } else group[target_group]->fill_queue.push_front(now_active); //push the active segment into FIFO queue
            ssd->TOTAL_GNUM[target_group] += 1; //increase the number of segments counter
            //TODO seg stamp changed
            ssd->seg_stamp[now_active] = stats->cur_wp; //not used yet
            ssd->fill_info[now_active] = 1; //mark as filled

            now_active = queue_pop_back(&fbqueue); //assign new active segment from free segment queue
            ssd->gnum_info[now_active] = target_group;
            ssd->active[target_group][0] = now_active; //active segment setting 
            ssd->active[target_group][1] = 0;
            ssd->seg_stamp[now_active]=-1;
        }

        if(ssd->itable[victim_idx+i] == false && ssd->oob[victim_idx+i].lba != -1){ //if the page in the victim segment is valid,
            local_copy += 1;
            copy_lba = ssd->oob[victim_idx+i].lba; //get lba of valid page
            if (victim_gnum == 0 || victim_gnum == 1) {
                hf_generate (ssd, copy_lba, victim_seg, group, hf_gen, 0);
            }
            new_ppa = now_active*ssd_spec->PPS + ssd->active[target_group][1]; //assign new ppa for valid page
            ssd->mtable[copy_lba] = new_ppa; //update mapping entry for valid page
            ssd->oob[new_ppa].lba = copy_lba; //record lba in the oob space 
            //setting segment age when the first page is written
            ssd->active[target_group][1] += 1; //update appending point of target group
        }
    }

    //update valid ratio and erase count information 
    //* add to ssd->VPG and ssd->EPG
    //* calculate segment age and add to HotFilter
    update_gc_results(ssd, stats, local_copy, victim_gnum, victim_seg);
    //update GC information
    stat_update(stats, local_copy, 1);

    //initialize segment information
    initialize_segment(ssd, victim_idx, victim_seg, victim_gnum);
    fbqueue.push_front(victim_seg);  //push the erased segment in the free segment queue

    return target_group;
}


int error_cmp(struct SSD *ssd, class GROUP *group[]){
    //return -1;
    printf("\n**Error comparing function**\n");
    double opt_vr[10];
    unsigned int new_group[10];

    int cur_ver = ssd->config_version;
    memcpy(opt_vr, ginfo->commit_vr, sizeof(double)*10);
    int loop_num = ssd_spec->GROUPNUM;
    if (ssd_spec->mida_on == 1) loop_num = ssd_spec->naive_start;
    //checking hot filter traffic, if lower than 0.75% of the first traffic, then turn off the hotfilter

    for(int i = 1; i < loop_num; i++){
        double real_vr = (double) ssd->err_VPG[i]/(double) ssd->err_EPG[i];
        printf("[group %d]'s real vs exp valid ratio:[%.4f vs %.4f]\n", i, real_vr, opt_vr[i]);
        if(opt_vr[i]-real_vr <= -0.1){ //error rate is more than 10%
        //if(((opt_vr[i]-real_vr)>=0.1)||(opt_vr[i]-real_vr <= -0.1)){ //error rate is more than 10%
            printf("**Detecting <Unknown Erea>**\n");
            return i;
        }
    }
    return -1;
}

//Group Split function
//old_last -> new_last
//!!!!! can't use this function when naive MiDA is on !!!!!
int group_split(struct SSD *ssd, struct STATS *stats, class GROUP *group[], int new_last){
    //there is no setup for naive MiDA mode
    if (ssd_spec->mida_on==1) {
        abort();
    }

    int old_last = ssd_spec->naive_start; //naive mida starting point
    printf("[ALERT!] GROUP SPLIT START <%d -> %d>!!\n", old_last, new_last);
    printf("********************************************\n");
    if (old_last == new_last) {
        abort();
    }
    int cnt = 0;

    //Group Structure Swapping
    GROUP* tmp = group[old_last];
    group[new_last] = tmp;
    ssd->TOTAL_GNUM[new_last] = ssd->TOTAL_GNUM[old_last];
    ssd->TOTAL_GNUM[old_last]=0;

    //making new group classes for new groups
    for (int i=old_last; i<new_last; i++) {
        group[i] = new GROUP;
        group[i]->size = 0;
        ssd->TOTAL_GNUM[i] = 0;
    }

    //Preparation of moving information of active segments of old last group
    if (ssd->active[new_last][0] != -1) {
        abort();
    }
    ssd->active[new_last][0] = ssd->active[old_last][0];
    ssd->active[new_last][1] = ssd->active[old_last][1];
    ssd->gnum_info[ssd->active[new_last][0]] = new_last;
    ssd->active[old_last][0]=-1;
    ssd->active[old_last][1]=0;

    //change the gnum_info value
    int old_total=0;
    int fill_size = group[new_last]->fill_queue.size();
    int mig_idx=0;
    int err_list[20] = {0,};
    int err_idx=0;
    for (int i=0;i<fill_size;i++) {
        mig_idx = group[new_last]->fill_queue.front();
        group[new_last]->fill_queue.pop_front();
        if (ssd->gnum_info[mig_idx] == old_last) {
            ssd->gnum_info[mig_idx] = new_last;
            old_total++;
        } else if (err_idx < 20){
            err_list[err_idx] = ssd->gnum_info[mig_idx];
            err_idx++;
        }
        group[new_last]->fill_queue.push_back(mig_idx);
    }
    if (old_total != fill_size) {
        for (int i=0;i<20;i++) printf("%d ", err_list[i]);
        printf("\n");
        abort();
    }
    if (ssd->TOTAL_GNUM[new_last] != old_total) {
        ssd->TOTAL_GNUM[new_last] = old_total;
        //TODO Check this!!!
        //abort();
    }

    //change the values
    ssd_spec->naive_start = new_last; //increase the starting point of naive mida
    ssd_spec->GROUPNUM = new_last+1; //increase of the number of groups

    //setting free active segment
    int t = 0;
    int ret_flag = 0;
    int ret;
    for (int i=old_last;i<new_last;i++) {
        if (ssd->active[i][0] != -1) {
            abort();
        }
        ssd->active[i][0]=-1;
        ssd->active[i][1]=0;
    }
    //if the old naive_start is 1, then need to make active segment to group 1
    if (old_last == 1) {
        ssd->active[1][0] = queue_pop_back(&fbqueue);
        ssd->gnum_info[ssd->active[1][0]] = 1;
    }

    printf("[ALERT!] GROUP SPLIT END!!\n");
    printf("********************************************\n");
    return 0;
}

//Active segment allocation for turning ON the naive mida
void active_split(struct SSD *ssd, struct STATS *stats, class GROUP *group[]){
    int gc_group;
    int ret;
    int gc_cnt = 0;
    while(fbqueue.size() < ssd_spec->FREENUM+(ssd_spec->NAIVE_N-1)){ //ensure the # of free segment for allocating the active segment for group of naive mida
        gc_group = do_gc_for_split(ssd, stats, group, ssd_spec->naive_start);
        if (gc_cnt == ILOOP) {
            merge_group(ssd,stats,group); //if there is no capability to provide free space, merge group
            if (ssd_spec->mida_on == 0) ret = mida_on_off(ssd, stats, group, ON); //If naive mida is OFF, Turn on the naive mida by unknown area
            gc_cnt = 0;
        }
        gc_cnt ++;
    }
    int new_start = ssd_spec->naive_start+1;
    //active segment allocation
    for (int i = 0; i < (ssd_spec->NAIVE_N-1); i++){
        ssd->active[new_start+i][0] = queue_pop_back(&fbqueue);
        ssd->active[new_start+i][1] = 0;
        ssd->gnum_info[ssd->active[new_start+i][0]] = new_start+i;
    }
    group[ssd_spec->naive_start]->size -= (ssd_spec->NAIVE_N-1); //Decrease the size of last group by active segments (newly allocated 3 segments)
    ret = gsize_check(ssd, group); //verify the size of each group
}

//Naive MiDA OFF
int mida_off(struct SSD *ssd, struct STATS *stats, class GROUP *group[]){

    int start_g = ssd_spec->naive_start;
    //merge the active segment of group 2 ~ 5 to group 1 queue
    for(int i = 0; i < (ssd_spec->NAIVE_N-1); i++){
        //start_g = naive_start (1)
        //+1 = active 1(group 1)
        int vic_active = ssd->active[start_g+1+i][0];
        if (vic_active == -1) continue;
        group[start_g]->fill_queue.push_front(vic_active);
        ssd->TOTAL_GNUM[start_g] += 1;
        ssd->seg_stamp[vic_active] = stats->cur_wp;
        ssd->fill_info[vic_active] = 1;
        //0 is HOT, naive mida groups are changing to group 1
        ssd->gnum_info[vic_active] = start_g;
        ssd->irtable[vic_active]=ssd_spec->PPS/2;
        ssd->active[start_g+1+i][0]=-1;
        ssd->active[start_g+1+i][1]=0;
    }
    return 0;
}

//Active Segment Merge (by Group merge)
int active_merge(struct SSD *ssd, struct STATS *stats, class GROUP *group[], int target_g, int victim_g){
    int copy_lba, new_ppa;
    int fill_idx = -1;
    int old_active = ssd->active[victim_g][0];
    if (old_active == -1) {
        return -1;
    }
    int victim_idx = old_active*ssd_spec->PPS;
    int old_active_fill = ssd->active[victim_g][1];
    int now_active = ssd->active[target_g][0];
    for (int i = 0; i < old_active_fill; i++){
        if(ssd->active[target_g][1] == ssd_spec->PPS){ // if active segment is full
            group[target_g]->fill_queue.push_front(now_active);
            ssd->TOTAL_GNUM[target_g] += 1; //increase the number of segments counter
            //TODO seg stamp
            ssd->seg_stamp[now_active] = stats->cur_wp; //not used yet
            ssd->fill_info[now_active] = 1; //mark as filled
            fill_idx = now_active;

            now_active = queue_pop_back(&fbqueue); //assign new active segment from free segment queue
            //TODO target_g?? victim_g???
            ssd->gnum_info[now_active] = target_g;
            ssd->active[target_g][0] = now_active; //active segment setting
            ssd->active[target_g][1] = 0;
            //ssd->seg_stamp[now_active]=stats->cur_wp;
        }

        if(ssd->itable[victim_idx+i] == false){ //if the page in the victim segment is valid,
            //stats->tmp_copy += 1; //increase copy counter
            copy_lba = ssd->oob[victim_idx+i].lba; //get lba of valid page
            new_ppa = now_active*ssd_spec->PPS + ssd->active[target_g][1]; //assign new ppa for valid page
            ssd->mtable[copy_lba] = new_ppa; //update mapping entry for valid page
            ssd->oob[new_ppa].lba = copy_lba; //record lba in the oob space
            ssd->active[target_g][1] += 1; //update appending point of target group
        }
    }
    for (int i = 0; i < ssd_spec->PPS; i++) { //erase segment
        ssd->itable[victim_idx+i] = 0; //invalid T/F clear
        ssd->oob[victim_idx+i].lba = -1; //oob clear
    }
    ssd->irtable[old_active] = 0; //invalid page counter clear
    ssd->gnum_info[old_active] = ssd_spec->MAXGNUM - 1; //group number clear, not used
    ssd->fill_info[old_active] = 0; //filled T/F clear
    ssd->seg_stamp[old_active] = -1;

    fbqueue.push_front(old_active);  //push the erased segment in the free segment queue
    ssd->active[victim_g][0]=-1;
    ssd->active[victim_g][1]=0;
    return fill_idx;
}

void hot_merge(struct SSD* ssd, struct STATS* stats, class GROUP *group[]) {
    int HOT=0;
    int G1=1;
    int move_seg;
    int hot_size = (int)(group[HOT]->fill_queue.size());
    for (int i=0;i<hot_size;i++) {
        move_seg = queue_pop_back(&group[HOT]->fill_queue);
        group[G1]->fill_queue.push_back(move_seg);
    }
    group[HOT]->size=0;
    return;
}

//Group Merge
//target_g is always naive_start (can't used in naive_mida groups)
void merge_group(struct SSD *ssd, struct STATS *stats, class GROUP *group[]){
    int gnum, trash, mig_idx, free_idx, move_idx, new_free_idx, old_free_idx, fill_size, victim_g;
    int cnt = 0;
    int mida_on = 0;
    int ret;
    int target_g = ssd_spec->naive_start;
    victim_g = target_g-1;

    if (target_g == 1) {
        abort();
    }

    printf("********************************************\n");
    printf("[ALERT!] GROUP MERGE START!! <%d-%d>\n", victim_g, target_g);
    fill_size = group[victim_g]->fill_queue.size();
    //TOTAL GNUM size check
    if (fill_size != ssd->TOTAL_GNUM[victim_g]) {
        abort();
    }
    //move segment to new group queue
    for (int i = 0; i < fill_size; i++){
        mig_idx = group[victim_g]->fill_queue.front();
        group[victim_g]->fill_queue.pop_front();
        ssd->gnum_info[mig_idx] += 1;
        group[target_g]->fill_queue.push_back(mig_idx);
        ssd->TOTAL_GNUM[target_g] += 1;
    }

    ret = active_merge(ssd, stats, group, target_g, victim_g);
    group[target_g]->size += group[victim_g]->size + 1; //+1: active segment

    //change the gnum info
    for(int i = 0; i < group[target_g]->fill_queue.size(); i++){
        int tmp_idx = queue_pop_back(&group[target_g]->fill_queue);
        ssd->gnum_info[tmp_idx] -= 1;
        group[target_g]->fill_queue.push_front(tmp_idx);
    }

    //move the active segment information
    for(int i = 0; i < ssd_spec->mida_on*(ssd_spec->NAIVE_N)+1; i++) {
        int tmp_gnum = ssd_spec->naive_start+i;
        int vic_gnum = tmp_gnum-1;
        int tmp_total = ssd->TOTAL_GNUM[tmp_gnum];
        int tmp_active_idx = ssd->active[tmp_gnum][0];
        int tmp_active_offset = ssd->active[tmp_gnum][1];

        ssd->TOTAL_GNUM[vic_gnum] = tmp_total;
        ssd->active[vic_gnum][0] = tmp_active_idx;
        ssd->active[vic_gnum][1] = tmp_active_offset;
        ssd->gnum_info[ssd->active[vic_gnum][0]] = vic_gnum;
        ssd->active[tmp_gnum][0] = -1;
        ssd->active[tmp_gnum][1] = 0;
        ssd->TOTAL_GNUM[tmp_gnum] = 0;
    }

    delete group[victim_g];
    group[victim_g] = group[target_g];
    group[target_g] = NULL;
    ssd_spec->GROUPNUM -= 1;
    ssd_spec->naive_start -= 1;
    printf("[ALERT!] GROUP MERGE END!! GROUPNUM: %d->%d\n", ssd_spec->GROUPNUM+1, ssd_spec->GROUPNUM);
    printf("GROUPNUM: %d\n", ssd_spec->GROUPNUM);
    printf("********************************************\n");

    //move hot segment to G1 queue (All groups are merged, so dont need to fix the HOT group size)
    if (ssd_spec->naive_start == 1) {
        hot_merge(ssd, stats, group);
        stats->err_flag=0;
        stats->err_write=0;
        ssd->config_version=0;
    }
    return;
}

int gsize_check(struct SSD* ssd, class GROUP *group[]){
    return 0;
    double gsize[20] = {0,};
    for (int i = 0; i < ssd_spec->SEGNUM; i++){
        if (ssd->gnum_info[i] < ssd_spec->GROUPNUM) {
            gsize[ssd->gnum_info[i]]++;
        }
    }
    for (int i=0;i<ssd_spec->GROUPNUM;i++) {
        if (ssd->active[i][0] != -1) gsize[i]--;
    }
    for (int i = 0; i < ssd_spec->MAXGNUM; i++){
        if (i < ssd_spec->GROUPNUM) {
            if (gsize[i] == 0) continue;
            if (ssd->TOTAL_GNUM[i] != (gsize[i])) {
                abort();
            }
            if ((i < ssd_spec->naive_start) && ((gsize[i]) != group[i]->fill_queue.size())) {
                abort();
            }
            ssd->TOTAL_GNUM[i] = gsize[i]-1.;
        } else ssd->TOTAL_GNUM[i] = 0;
    }
    return 0;
}

int move_margin_seg(struct SSD* ssd, struct STATS *stats, class GROUP *group[], int victim_g, int target_g, int size){
    int move_idx;
    int prev_group = 0;
    for (int i = 0 ; i < size; i++){
        move_idx = queue_pop_back(&group[victim_g]->fill_queue);
        prev_group = ssd->gnum_info[move_idx];
        ssd->gnum_info[move_idx] = target_g;
        group[target_g]->fill_queue.push_front(move_idx);

        ssd->TOTAL_GNUM[target_g]++;
        ssd->TOTAL_GNUM[victim_g]--;

    }

    return 0;
}

int change_config_2(struct SSD* ssd, struct STATS *stats, class GROUP *group[]){
    int old_gnum = ssd_spec->naive_start+1;
    int new_gnum = ginfo->gnum;
    ginfo->M_gnum = ginfo->gnum;
    int margin_fseg = 0;
    //prev group number is bigger than the new one
    if (new_gnum < old_gnum) {
        for (int i = 0; i < old_gnum-new_gnum; i++){
            merge_group(ssd, stats, group);
        }
        margin_fseg = fbqueue.size()-ssd_spec->FREENUM;
    } else if (new_gnum > old_gnum){
        group_split(ssd, stats, group, new_gnum-1);
    }

    for (int i=0;i<new_gnum;i++) {
        if (group[i]== NULL) {
            abort();
        }
        group[i]->size = ginfo->gsize[i]-1;
        if (i==1) group[i]->size -= 1; //active segment
        if (i != ssd_spec->naive_start && group[i]->size < ssd->TOTAL_GNUM[i]) {
            move_margin_seg(ssd, stats, group, i, ssd_spec->naive_start, ssd->TOTAL_GNUM[i]-group[i]->size);
        }
    }

    return margin_fseg;
}

//checking if ginfo is valid
//if ginfo is valid, then compare the WAF and apply the new configuration
int checking_ginfo(struct SSD* ssd, struct STATS* stats, class GROUP *group[]) {
    if ((ginfo->valid == false)) return 0;
    int ret=0;
    int APPLY_FLAG=0;
    //check if the new WAF is lower than the current WAF
    APPLY_FLAG = check_applying_config(ssd, stats, group);
    if (APPLY_FLAG==1) {
        //if the current model has a lot of err, then double the next model's time window size
        if (mmodel->tmp_err_cnt != 0) mmodel->err_cnt++;
        else mmodel->err_cnt=0;
        if (mmodel->err_cnt >= 2) {
            model_resizing_timewindow(mmodel, mmodel->time_window*2);
            ginfo->err_stat=true;
        } else ginfo->err_stat=false;
    } else if (mmodel->tmp_err_cnt != 0) mmodel->err_cnt++;
    mmodel->tmp_err_cnt=0;
    //turn on the new model
    if (APPLY_FLAG==1) {
        //naive mida off
        if (ssd_spec->mida_on==1) ret=mida_on_off(ssd, stats, group, OFF);
        ginfo->commit_g=0;
        //merge or split
        ret=change_config_2(ssd, stats, group);
        ret=gsize_check(ssd, group);
        memcpy(ginfo->commit_vr, ginfo->vr, sizeof(double)*10);
        ssd->config_version=1;
        ssd->gnum = ginfo->gnum;
        stats->commit_bit=0;

        //reset the err flag
        stats->err_flag=0;
        stats->err_write=0;
    }
    model_initialize(mmodel, true);
    ginfo->valid=false;
    APPLY_FLAG=0;
    return 1;
}

//check if err handling is needed
void checking_err_time(struct SSD* ssd, struct STATS* stats, class GROUP *group[]) {
    if (stats->err_flag==1) {
        stats->err_write += 1;
    }
    if (stats->err_write > (long)mmodel->real_tw*0.3) {
        err_handling(ssd, stats, group);
        stats->err_write=0;
    }
    return;
}

void modeling_check(struct SSD* ssd, struct STATS* stats, class GROUP *group[], int user_group, int lba) {
    //TODO change the condition to checking interval
    if (mmodel->model_on) {
        int res = check_time_window(mmodel, WRITE);
        if (res) check_interval(mmodel, lba, WRITE);
    }
    //change the group configuration when ginfo is valid
    checking_ginfo(ssd, stats, group);

    //checking err_write, and handling err
    checking_err_time(ssd, stats, group);

    return;
}

void checking_config_apply(struct SSD* ssd, struct STATS* stats, class GROUP *group[], int gc_group) {
    if (ginfo->commit_g == 19) return;
    for (int i=0;i<ssd_spec->naive_start;i++) {
        if (stats->commit_bit && (int)pow(2, i)) ginfo->commit_g=i;
        else break;
    }
    if (stats->commit_bit == pow(2,ssd_spec->naive_start+1)-1) {
        stats->err_flag=1;
        ginfo->commit_g=19;
    }
    return;
}

int infinite_gc_handling(struct SSD* ssd, struct STATS* stats, class GROUP *group[], int gc_group) {
    //add model error count
    if (ssd->config_version==1 && gc_group <= (ssd->gnum/2)) {
        mmodel->tmp_err_cnt++;
    }
    if(ginfo->commit_g != 19) {
        gc_group = ginfo->commit_g;
    }
    //merge_group
    if (ssd_spec->naive_start != 1){
        if (gc_group >= ssd_spec->naive_start) gc_group = ssd_spec->naive_start;
        gc_group -= 1;
        for (int i=ssd_spec->naive_start;i!=gc_group;i--) {
            if (ssd_spec->naive_start == 1) break;
            merge_group(ssd, stats, group);
        }
        if (ginfo->commit_g != 19 && ginfo->commit_g >= ssd_spec->naive_start) {
            stats->err_flag=1;
            ginfo->commit_g=19;
        }

        if (ssd_spec->mida_on==0) {
            mida_on_off(ssd, stats, group, ON);
        }
        return 1;
    } else {
        ssd->config_version = 0;
        return 0;
    }
}

//just configure WAF
int check_applying_config(struct SSD *ssd, struct STATS *stats, class GROUP* group[]){
    double new_waf = ginfo->WAF;
    double cur_waf = stats->tmp_waf;

    stats->err_write = 0;
    stats->err_flag = 0;
    for(int i = 0; i < ssd_spec->MAXGNUM; i++){
            ssd->err_VPG[i] = 0.;
            ssd->err_EPG[i] = 0.;
    }

    if (cur_waf < new_waf) {
//      printf("NOT NEED TO CHANGE:: NEW: %.3f, CUR: %.3f\n", new_waf, cur_waf);
        return 0;
    }
    else {
//      printf("NEED TO CHANGE:: NEW: %.3f, CUR: %.3f\n", new_waf, cur_waf);
        return 1;
    }
}

int mida_on_off(struct SSD *ssd, struct STATS *stats, class GROUP* group[], int SWITCH){
    int ret;
    if (SWITCH == ssd_spec->mida_on) return 0;
    if (SWITCH){ // MiDA ON
        ssd_spec->mida_on = 1;
        ssd_spec->GROUPNUM = ssd_spec->naive_start+ssd_spec->NAIVE_N;
    }
    else { // MiDA OFF
        int cur_gnum = ssd_spec->naive_start;
        int seg_idx;
        int tmp_size;
        //collect naive Mida's active segments (just add to the fill_queue)
        ret = mida_off(ssd, stats, group);

        ssd_spec->mida_on = 0;

        //change group number information (naive_start+1, 2, 3, 4 and 15 -> naive_start)
        //TODO group 15 segment in naive_queue is converted to naive_start group. is that OK?
        seg_idx=0;
        tmp_size = group[ssd_spec->naive_start]->fill_queue.size();
        for (int i=0;i<tmp_size;i++) {
            seg_idx = group[ssd_spec->naive_start]->fill_queue.front();
            group[ssd_spec->naive_start]->fill_queue.pop_front();
            if (ssd->gnum_info[seg_idx] > ssd_spec->naive_start) {
                ssd->gnum_info[seg_idx] = ssd_spec->naive_start;
            }
            group[ssd_spec->naive_start]->fill_queue.push_back(seg_idx);
        }
        for (int i=1;i<ssd_spec->NAIVE_N;i++) {
            ssd->TOTAL_GNUM[ssd_spec->naive_start] += ssd->TOTAL_GNUM[ssd_spec->naive_start+i];
            ssd->TOTAL_GNUM[ssd_spec->naive_start+i]=0;
        }

        //move group 0(HOT) segments to group[0]
        seg_idx=0;
        tmp_size = group[ssd_spec->naive_start]->fill_queue.size();
        int g0_cnt=0;
        if (ssd_spec->naive_start == 1) {
            for (int i=0;i<tmp_size;i++) {
                seg_idx = group[1]->fill_queue.front();
                group[1]->fill_queue.pop_front();
                if (ssd->gnum_info[seg_idx] == 0) {
                    group[0]->fill_queue.push_back(seg_idx);
                    g0_cnt++;
                } else group[1]->fill_queue.push_back(seg_idx);
            }
            if (ssd->TOTAL_GNUM[0] != g0_cnt) {
                abort();
            }
        }


        ssd_spec->GROUPNUM = ssd_spec->naive_start+1;
        ret = gsize_check(ssd, group);
    }
    return 0;
}

void err_handling(struct SSD* ssd, struct STATS* stats, class GROUP* group[]){
    stats->err_flag = 1;
    int unknown_num = -1;
    int diff = 0;
    unknown_num = error_cmp(ssd, group);
    stats->unknown_cnt = 0;
    if (unknown_num != -1){
        printf("%dth group is Unknown Area\n", unknown_num);
        diff = ssd_spec->naive_start - unknown_num;
        if (!ssd_spec->mida_on) mida_on_off(ssd, stats, group, ON);
        for (int i = 0; i < diff; i++){
            merge_group(ssd, stats, group);
        }
        if ((unknown_num <= (ginfo->M_gnum/2))) {
            mmodel->tmp_err_cnt++;
        }
    }
    for (int i = 0; i < ssd_spec->MAXGNUM; i++){
        ssd->err_VPG[i] = 0.;
        ssd->err_EPG[i] = 0.;
    }
}


