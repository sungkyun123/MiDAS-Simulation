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
#include <unistd.h>

#include "model.h"
#include "ssd_config.h"
#include "queue.h"
#include "hf.h"
#include "algorithm.h"
#include "ssdsimul.h"

using namespace std;

//progress when apply group size

#define OFF 0 // naive_mida OFF
#define ON 1 // naive_mida ON

#define ILOOP 25 //When the # of recursive GC is over ILOOP, merge the last two groups

#define GB_P (1024*1024/4)
#define SW_GENERATE 1 // 0: No 1: Yes

uint32_t utilization=0;
struct SSD_SPEC *ssd_spec;
struct SSD *ssd;
struct STATS *stats;

char *workload;
char *c_workload;
char *vs_policy;

deque<int> fbqueue; //global free block management queue

mini_model *mmodel;
G_INFO *ginfo;
MODEL_Q *model_q;

struct HotFilter *hf_gen;

double age_function(struct SSD *ssd, struct STATS *stats, int seg_idx){
	if (ssd->seg_stamp[seg_idx] == -1) {
		abort();
	}	
	long age_real = stats->cur_wp - ssd->seg_stamp[seg_idx];
	return sqrt(age_real);
}

int CB_VS(struct SSD *ssd, struct STATS *stats, class GROUP *group[], int gc_group){
	int cand_segment = -1;
	double cand_value = -1;
	int comp_segment = 0;
	double comp_value = -1.;
	double age_value = 0;
	int queue_size = (int)(group[gc_group]->fill_queue.size());
	if ((gc_group < ssd_spec->naive_start) && (queue_size != ssd->TOTAL_GNUM[gc_group])) {
		abort();
	}
	for(int i = 0; i < queue_size; i++){ //pop segment from the FIFO queue and compared with cand semgnet
		comp_segment = queue_pop_back(&group[gc_group]->fill_queue);
		age_value = age_function(ssd, stats, comp_segment);
		if (ssd->fill_info[comp_segment]==0) {
		       abort();
		}	       
		if (ssd->irtable[comp_segment] == ssd_spec->PPS) {
			if (cand_segment != -1) group[gc_group]->fill_queue.push_front(cand_segment);
			return comp_segment;
		}
		//age_value = 1;
		comp_value = (double)(ssd->irtable[comp_segment])*(double)age_value/((double)ssd_spec->PPS - (double)ssd->irtable[comp_segment]);
		if (comp_value > cand_value) { //if popped one has higher invalid ratio than the candidate segment one;
			if (cand_segment != -1) group[gc_group]->fill_queue.push_front(cand_segment); //enqueue the candidate one into FIFO queue
			cand_segment = comp_segment; //update the candidate to compared one
			cand_value = comp_value;
		}
		else{
			group[gc_group]->fill_queue.push_front(comp_segment); //else, enqueue the compare one into FIFO queue again
		}
	}
		return cand_segment; //select the vctim segment with highest invalid ratio in the GC target group
}


//FIFO policy
int FIFO_VS(class GROUP *group[], int gc_group){
	return queue_pop_back(&group[gc_group]->fill_queue); //pop the segment from the victim group FIFO queue
}

//Greedy policy
int GREEDY_VS(struct SSD *ssd, class GROUP *group[], int gc_group){ 
	int cand_segment = queue_pop_back(&group[gc_group]->fill_queue); //initialize candidate victim segment 
	int cand_ir = ssd->irtable[cand_segment]; //identify the invalid ratio of candidate victim segment
	int comp_segment; //for compare segment
	int comp_ir; //for compare segment 
	int queue_size = (int)(group[gc_group]->fill_queue.size());
	if ((gc_group < ssd_spec->naive_start) && (queue_size != ssd->TOTAL_GNUM[gc_group])) {
		abort();
	}
	for(int i = 0; i < queue_size; i++){ //pop segment from the FIFO queue and compared with cand semgnet
		comp_segment = queue_pop_back(&group[gc_group]->fill_queue);
		comp_ir = ssd->irtable[comp_segment];
		if (cand_ir < comp_ir) { //if popped one has higher invalid ratio than the candidate segment one;
			group[gc_group]->fill_queue.push_front(cand_segment); //enqueue the candidate one into FIFO queue
			cand_segment = comp_segment; //update the candidate to compared one
			cand_ir = comp_ir;
		}
		else{
			group[gc_group]->fill_queue.push_front(comp_segment); //else, enqueue the compare one into FIFO queue again
		}
	}
	return cand_segment; //select the vctim segment with highest invalid ratio in the GC target group
}

//Victim selection interface
int victim_selection(struct SSD *ssd, struct STATS *stats, class GROUP *group[], int gc_group){
	int victim_seg;
	if (gc_group < ssd_spec->naive_start) {
		victim_seg = FIFO_VS(group, gc_group);
	} else {
		if (ssd_spec->vs_policy == 0) victim_seg = FIFO_VS(group, gc_group);
		else if (ssd_spec->vs_policy == 1) victim_seg = GREEDY_VS(ssd, group, gc_group);
		else if (ssd_spec->vs_policy == 2) {
			victim_seg = CB_VS(ssd, stats, group, gc_group);
		}
	}
	return victim_seg;
}


// update the results of GC
// e.g., valid ratio, erase count, segment age in HotFilter
void update_gc_results(struct SSD* ssd, struct STATS* stats, int local_copy, int victim_gnum, int victim_seg) {
	//update segment age when group number is 0
	if ((hf_gen->make_flag==1) && (victim_gnum==0)) {
		double seg_age = (double)(stats->cur_wp - ssd->seg_stamp[victim_seg]);
		hf_gen->seg_age += seg_age;
		hf_gen->seg_num++;
		hf_gen->G0_vr_sum += (double)(local_copy)/(double)(ssd_spec->PPS);
		hf_gen->G0_vr_num++;
	}
	
	if (stats->progress >= 10){
		ssd->VPG[victim_gnum] += (double)(local_copy)/(double)(ssd_spec->PPS); //calculate valid ratio per group		    
		ssd->EPG[victim_gnum] += 1; //increase erase count of group
	}
	//add GC information for error check
	int cur_ver = ssd->config_version;
	if (cur_ver >= 0) {
		if (ginfo->commit_g == 19 && (stats->err_flag == 1)) {
			//stats->err_flag = 1;
			ssd->err_VPG[victim_gnum] += (double)(local_copy)/(double)(ssd_spec->PPS); //calculate valid ratio per group	    
			ssd->err_EPG[victim_gnum] += 1; //increase erase count of group
		}
	}
	ssd->TMP_VPG[victim_gnum] += (double)(local_copy)/(double)(ssd_spec->PPS);
	ssd->TMP_EPG[victim_gnum] += 1;
	ssd->DEBUG_VPG[victim_gnum] += (double)(local_copy)/(double)(ssd_spec->PPS);
	ssd->DEBUG_EPG[victim_gnum] += 1;
	
	return;
}

//initialize segment information when the segment is erased
void initialize_segment(struct SSD *ssd, int victim_idx, int victim_seg, int victim_gnum) {
	for (int i = 0; i < ssd_spec->PPS; i++) { //erase segment
		ssd->itable[victim_idx+i] = 0; //invalid T/F clear
		ssd->oob[victim_idx+i].lba = -1; //oob clear
	}
	ssd->irtable[victim_seg] = 0; //invalid page counter clear
	ssd->gnum_info[victim_seg] = ssd_spec->MAXGNUM - 1; //group number clear, not used
	ssd->fill_info[victim_seg] = 0; //filled T/F clear
	ssd->TOTAL_GNUM[victim_gnum] -= 1; //decrease the number of segment in the group 
	ssd->seg_stamp[victim_seg] = -1;
	return;
}

int do_gc(struct SSD *ssd, struct STATS *stats, class GROUP *group[], int gc_group){
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

//find the gc victim group based on the designated group size, 
int gc_victim_group(struct SSD *ssd, class GROUP *group[]) {
	//the case that naive_mida mode only
	if (ssd_spec->naive_start == 1) return ssd_spec->naive_start;
	
	//HOT group size is fixed
	for (int i=0;i<ssd_spec->naive_start;i++) {
		if (group[i]->size < ssd->TOTAL_GNUM[i]) {
			return i;
		}
	}
	return ssd_spec->naive_start;
}

//do GC while the fbqueue size over the FREENUM size
//but return 1 when the infinite GC occurs.
int GC(struct SSD* ssd, struct STATS *stats, class GROUP *group[]) {
	int t=0;
	int gc_group=0;
	while(fbqueue.size() < ssd_spec->FREENUM) {
		t++;
		gc_group = gc_victim_group(ssd, group);
		stats->commit_bit = stats->commit_bit | (int)pow(2, gc_group);
		gc_group = do_gc(ssd, stats, group, gc_group);
		if (t % ILOOP == 0 && t>1) {
			return gc_group;
		}
	}
	checking_config_apply(ssd, stats, group, gc_group);
	return 0;
}

int write(int lba, struct SSD *ssd, struct STATS *stats, class GROUP *group[], bool sw){
	int user_group = 0;

	//TODO hf_generate timing
	//if hf_generate first, then the lba will be assigned to HOT group immediately when the past lba has the short interval
	
	int old_ppa = ssd->mtable[lba]; //get the previous ppa of lba
	int old_seg = old_ppa/ssd_spec->PPS; //get the previous segment of lba
	int old_group = ssd->gnum_info[old_seg];
	if(old_ppa != -1){ //if there is old data of lba in SSD:
		ssd->itable[old_ppa] = 1; //mark as invalid
		ssd->irtable[old_seg] ++; //increase invalid count of segment
		hf_generate(ssd, lba, old_seg, group, hf_gen, 1);
	} else utilization++;
	

	user_group = hf_check(lba, hf_gen);
	int active_old = ssd->active[user_group][0]; //index of now active block
	if (active_old == -1) {
		abort();
	}
	int active_new; //index of new active block
	int local_copy = 0;
	int gc_group;
	int real_gc_group;
	int ret;

	if(ssd->active[user_group][1] == ssd_spec->PPS){ //if active block is full:
		//separate when the group is naive_mida mode or MiDASS mode
		if (user_group > ssd_spec->naive_start) {
			abort();
			group[ssd_spec->naive_start]->fill_queue.push_front(active_old);
		} else if (ssd_spec->naive_start == 1) group[ssd_spec->naive_start]->fill_queue.push_front(active_old);
		else group[user_group]->fill_queue.push_front(active_old);

		ssd->TOTAL_GNUM[user_group] += 1;
		//TODO seg stamp
		ssd->seg_stamp[active_old] = stats->cur_wp; //record the sealed time of active block
		ssd->fill_info[active_old] = 1; //set 1 by filled state of active block

		active_new = queue_pop_back(&fbqueue); //assign new active block to group from global free queue
		ssd->active[user_group][0] = active_new; // set new active block
		ssd->active[user_group][1] = 0;
		ssd->gnum_info[active_new] = user_group;
		//ssd->seg_stamp[active_new] = stats->cur_wp;
	}


	int new_ppa = ssd->active[user_group][0]*ssd_spec->PPS+ssd->active[user_group][1];
	
	if (ssd->mtable[lba] == -1) stats->vp += 1;
	ssd->mtable[lba] = new_ppa; //update mapping table of lba
	ssd->oob[new_ppa].lba = lba; //record lba and group number in oob space of this page
	ssd->oob[new_ppa].gnum = user_group;
	ssd->active[user_group][1] ++; //increase appending point of active block

	//handling modeling functions
	modeling_check(ssd, stats, group, user_group,lba);
	
	gc_group = user_group;
	
	//have to do GC!!!
	ret = GC(ssd, stats, group);
	if (ret != 0) {
		infinite_gc_handling(ssd, stats, group, ret);
		ssd_spec->INF_GC_CNT++;
	}else ssd_spec->INF_GC_CNT=0;

	//infinite GC occurred continuously
	if (ssd_spec->INF_GC_CNT>=100) {
		abort();
	}
	
	return 0;
}

int trim(int lba, struct SSD *ssd, struct STATS *stats){
	int trim_ppa = ssd->mtable[lba];

	if (ssd->mtable[lba] == -1) return 0;

	if (trim_ppa != -1) {
		utilization--;	
		stats->vp--;
		ssd->itable[trim_ppa] = 1;
		ssd->mtable[lba] = -1;
		ssd->irtable[(int)(trim_ppa/ssd_spec->PPS)] += 1;
	}
	if (mmodel->model_on) {
		int res = check_time_window(mmodel,REMOVE);
		if (res) check_interval(mmodel,lba, REMOVE);
	}
	return 0;
}
char workload_name[64]=""; // move to here for valid ratio graph by soyoung

//display the information of simulation for each 1% progress
int simul_info(struct SSD *ssd, struct STATS *stats, class GROUP* group[], int time_gap){

	double waf, tmp_waf;

	waf = (double)(stats->cur_wp+stats->copy)/(double)(stats->cur_wp);
	tmp_waf = (double)(stats->tmp_wp+stats->tmp_copy)/(double)(stats->tmp_wp);
	stats->tmp_waf = tmp_waf;
	stats->util = (double)stats->vp/(double)ssd_spec->LBANUM;
	if (stats->progress >= 20) {
		stats->fourty_waf += tmp_waf;
	}
	printf("\n[Progress: %dGB] WAF: %.3f, TMP_WAF: %.3f, Utilization: %.3f\n", (int)(time_gap)*stats->progress, waf, tmp_waf, stats->util);
	for (int i = 0; i < ssd_spec->GROUPNUM; i++){
		if(i == ssd_spec->naive_start){
			printf("--------------------------\n");
		}
		ssd->CUR_VPG[i] = ssd->TMP_VPG[i]/ssd->TMP_EPG[i];
		if (i==0) printf("HOT[%d]: %.4f (ERASE: %d) ", (int)(ssd->TOTAL_GNUM[i]), ssd->CUR_VPG[i], (int)ssd->TMP_EPG[i]);
		else printf("GROUP %d[%d]: %.4f (ERASE: %d) ", i, (int)(ssd->TOTAL_GNUM[i]), ssd->CUR_VPG[i], (int)ssd->TMP_EPG[i]);
		if (ssd_spec->naive_start != 1 && ssd_spec->naive_start >= i) {
			printf(" (Q size: %ld) (desig_size: %d)\n", group[i]->fill_queue.size(), group[i]->size);
		} else printf("\n");
	}
	printf("\n");
	int ret;
	stats->progress += 1;
	stats_clear(ssd, stats);
	sleep(5);
	return 0;
}

//request preprocess into [type, lba, io_size]
void req_processing(char* raw_req, struct REQUEST *req, struct STATS *stats){
	req->timestamp = stats->checktime + (double)(atof(strtok(raw_req, " \t")));
	req->type = (int)(atoi(strtok(NULL, " \t")));
	req->lba = (int)(atoi(strtok(NULL, " \t")));
	req->io_size = (int)(atoi(strtok(NULL, " \t")));
	req->stream=0;
	req->sw=false;
}

void req_processing_sw(struct REQUEST *req, struct STATS *stats, int lba){
	req->timestamp = 0.;
	req->type = 1;
	req->lba = lba;
	req->io_size = 4096;
	req->stream=0;
	req->sw=true;
}


//request process
int submit_io(struct SSD *ssd, struct STATS *stats, class GROUP* group[], struct REQUEST *req){
	int ret;
	int req_num = (int)(req->io_size/ssd_spec->PGSIZE); //# of splited request into page unit
	// type 0: read, ret->0
	// type 1: write, ret-># of write page
	// type 2: ?, ret->0
	// type 3: discard (trim), ret-> -(# of trim page)
	if(req->type == 0 || req->type == 2) ret = 0; //no write or trim command
	else if (req->type == 1) {
		for (int i = 0; i < req_num; i++){
			write(req->lba+i, ssd, stats, group, req->sw); //write for request io size
		}
		ret = req_num;
	}
	else if (req->type == 3) {
		for (int i = 0; i < req_num; i++){
			ret = trim(req->lba+i, ssd, stats); //discard for request io size
		}
		ret = 0;
	}
	return ret;
}

//Clear temp variable in stats when the time window reaches
int stats_clear(struct SSD* ssd, struct STATS* stats){
	stats->tmp_wp = 1;
	stats->tmp_copy = 0;
	for(int i = 0; i < ssd_spec->GROUPNUM+1; i++){
		ssd->TMP_VPG[i] = 0;
		ssd->CUR_VPG[i] = 0;
		ssd->TMP_EPG[i] = 0;
	}
	ssd->TMP_VPG[15]=0;
	ssd->CUR_VPG[15] = 0;
	ssd->TMP_EPG[15] = 0;
	ssd->TMP_VPG[16]=0;
	ssd->CUR_VPG[16] = 0;
	ssd->TMP_EPG[16] = 0;
	return 0;
}

bool run_signal = false;
//stat update function
void stat_update(struct STATS* stats, int ret, int type){
	if (type == 0){ //positive or 0: from request
		stats->cur_wp += ret;
		if (run_signal) stats->load_wp += ret;
		stats->tmp_wp += ret;
		stats->m_wp += ret;
		stats->cur_req++; //increase cur_req count by complated request 
	}//ret is positive, update for write request (ret: # of write page)
	else if (type == 1){ //negative: from GC
		stats->copy += ret;
		stats->tmp_copy += ret;
		stats->m_copy += ret;
		stats->erase ++; //GC
	}//else: update for GC request(ret: # of copy page)
	else if (type == 2) {
		stats->cur_req++;
	}
}

int ssd_simulation_sw(struct SSD* ssd, struct STATS* stats, class GROUP* group[], char *workload, int group_size[]){
	
	struct REQUEST *req;
	int ret;
	int type;
	req = (struct REQUEST*)malloc(sizeof(struct REQUEST));
	int time_gap = ssd_spec->dev_gb;
	for (int i = 0; i < ssd_spec->LBANUM; i++){
		if (stats->cur_wp%(time_gap*GB_P)==3) simul_info(ssd, stats, group, time_gap);
		req_processing_sw(req,stats,i);
		ret = submit_io(ssd, stats, group, req); //submit io request
		if (ret == 0) type = 2;
		else if (ret > 0) type=0;
		stat_update(stats, ret, type);
		if (stats->cur_wp % (100*GB_P) == 1 && stats->cur_req != 0) {
			stats->m_wp = 1;
			stats->m_copy = 0;
		}
	}
	return 0;
}

int ssd_simulation(struct SSD* ssd, struct STATS* stats, class GROUP* group[], char* workload, int group_size[]){

	char raw_req[128]; //for workload line of trace file
	char load_mark[128] = "LOAD_END\n"; //for workloads with LOAD
	struct REQUEST *req;
	int ret;

	req = (struct REQUEST*)malloc(sizeof(struct REQUEST)); //request infomation structure
	FILE *wktrace = fopen(workload, "r"); //open trace file
	int type = 0;
	uint32_t sw_slicer=0;
	int flag = 0;
	int sub_flag = 1;

	int time_gap=ssd_spec->dev_gb;

	while(fgets(raw_req, 128, wktrace)){
		if(!strcmp(raw_req, load_mark)) {
			stats->tot_req -= 1; 
			run_signal=true;
			continue;
		}//load part check
		//print simulation
		if (stats->cur_wp%(time_gap*GB_P)==3) simul_info(ssd, stats, group, time_gap);
		
		req_processing(raw_req, req, stats); // request processing ([type, lba, io_size, stream, timestamp])
		stats->cur_time = req->timestamp;
		ret = submit_io(ssd, stats, group, req); //submit io request
		if (ret == 0) type = 2;
		else if (ret > 0) type=0;
		stat_update(stats, ret, type);
		
		if (stats->cur_wp % (100*GB_P) == 0 && stats->cur_req != 0) {
			stats->m_wp = 1;
			stats->m_copy = 0;
			gsize_check(ssd, group);
		}
	}
	stats->checktime = stats->cur_time;
	return 0;	
}

int display_result(struct SSD *ssd, struct STATS *stats, char *workload, char *vs_policy){
	string technique;
	technique = "MiDAS";

	printf("Simulation setup\n============================================\n");
	printf("Technique: %s\n", technique.c_str());
	printf("workload: %s\n# of groups: %d\nvictim selection: %s\n", workload, ssd_spec->GROUPNUM, vs_policy);
	printf("write size: %fTB\t write after load size: %fTB\n", (double)(stats->cur_wp/262144)/1024.0, (double)(stats->load_wp/262144)/1024.0);
	printf("============================================\n");
	printf("\n Experimental result\n");
	printf("============================================\n");
	double waf = (double)(stats->cur_wp+stats->copy)/(double)(stats->cur_wp);
	double tmp_waf = (double)(stats->tmp_wp+stats->tmp_copy)/(double)(stats->tmp_wp);
	stats->fourty_waf += tmp_waf;
	printf("USERWRITE: %llu\n", stats->cur_wp);
	printf("COPYWRITE: %llu\n", stats->copy);
	printf("TOTAL WAF: %.3f\n", waf);
	for (int i = 0; i < ssd_spec->GROUPNUM; i++){
		if (i == 0) printf("HOT[%d]: %.4f (ERASE: %d)\n", (int)(ssd->TOTAL_GNUM[i]+1), ssd->VPG[i]/ssd->EPG[i], (int)ssd->EPG[i]);
		else printf("GROUP %d[%d]: %.4f (ERASE: %d)\n", i, (int)(ssd->TOTAL_GNUM[i]+1), ssd->VPG[i]/ssd->EPG[i], (int)ssd->EPG[i]);
	}
	return 0;
}

int display_debug(struct SSD *ssd, char *workload) {
	printf("\n========================================\n");
	printf("DEBUGGING INFO\n");
	for (int i = 0; i < ssd_spec->GROUPNUM; i++){
		printf("GROUP %d[%d]: %.4f (ERASE: %d) (tmp erase: %d)\n", i, (int)(ssd->TOTAL_GNUM[i]+1), ssd->DEBUG_VPG[i]/ssd->DEBUG_EPG[i], (int)ssd->DEBUG_EPG[i], (int)ssd->TMP_EPG[i]);
	}	
	printf("========================================\n");
	return 0;
}

void sig_handler(int signum) {
	display_result(ssd, stats, workload, vs_policy);
	display_debug(ssd, workload);
	exit(0);
}

int main(int argc, char** argv){
	if (argc < 2) {
		printf("./midas <Workload> <victim selection (fifo/greedy/cost-benefit)> <Logical device size (GB)> <Segment size (MB)>\n");
		abort();
	}
	setbuf(stdout, NULL);
	workload = argv[1]; //workload
	vs_policy = argv[2]; //victim selection policy
	int dev_gb = (int)(atoi(argv[3])); //GB
	int seg_mb = (int)(atoi(argv[4])); //MB
	int gnum = 1; //initial # of groups
	int ar2 = gnum;
	int tw = 4; //time window size (tw x device size)  
	double write_size = dev_gb*tw;
	int mida_on = 1; //1: run MiDA when starting
	int naive_start;
	if (mida_on == 1) {
		naive_start = gnum; //1
		gnum += 4;
	}
	else naive_start = gnum-1;
	printf("here\n");
	printf("======================SIMULATION SETUP=====================\n");
	printf("===========================================================\n");
	signal(SIGINT, sig_handler);
	setbuf(stdout, NULL);
	setbuf(stderr, NULL);

	printf("Workload: %s\n# of groups: %d\nvictim selection: %s\n", workload, gnum, vs_policy);
	int queue_gnum = naive_start+1; //2, # of Groups with Queue management 
	int group_size[queue_gnum]; 
	int sum;
	int pc_flag = 0;
	
	group_size[queue_gnum-1] = -1; //setup after ssd_init
	
	ssd= (struct SSD*)malloc(sizeof(struct SSD)); //ssd structure malloc
	stats = (struct STATS*)malloc(sizeof(struct STATS)); //data structure malloc

	GROUP *group[20];

	int intpolicy = policy_to_int(vs_policy);

	//simulation initialization	
	ssd_init(ssd, group, gnum, intpolicy, naive_start, dev_gb, seg_mb);
	group_init(ssd ,group, queue_gnum, mida_on, group_size);
	stats_init(stats);
	hf_init(ssd, &hf_gen);
	stats->tw = write_size;
	mmodel=(mini_model*)malloc(sizeof(mini_model));
	model_create(mmodel, write_size, 0);

	if (SW_GENERATE == 1) ssd_simulation_sw(ssd, stats, group, workload, group_size);
	hf_gen->make_flag=1;
	hf_gen->use_flag=1;
	ssd_simulation(ssd, stats, group, workload, group_size);

	display_result(ssd, stats, workload, vs_policy);
	pthread_join(mmodel->thread_id, NULL);
	model_destroy(mmodel);
}
