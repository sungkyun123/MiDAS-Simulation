#include "model.h"
#include "ssd_config.h"
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <boost/math/special_functions/lambert_w.hpp>
#include <cmath>

//extern algorithm page_ftl;
using boost::math::lambert_w0;
using boost::math::lambert_wm1;

extern uint32_t utilization;
extern struct SSD_SPEC *ssd_spec;
extern char workload_name[64];

extern struct HotFilter *hotfilter;
extern G_INFO *ginfo;
extern MODEL_Q *model_q;
extern struct SSD *ssd;
#define RESIZING 1 
#define DES 2//desnoyer, DES=2 means Greedy
#define ADD_TWO 0 //+=2
#define REAL_BENCH false

#define DES_FIX 1
#define HOT_POSITION 0

#define FIRST_NEW 0

/* calculate WAF value using valid ratio per groups
 * values = valid ratio per groups
 * return: calculated WAF
 */

double dqueue_pop_back(deque<double> *dqueue) {
	double ret = dqueue->back();
	dqueue->pop_back();
	return ret;
}

/*initialize model*/
void model_create(mini_model *mmodel, int write_size, int type) {
	/* sampling info, unit size, time window size */
	model_q = new MODEL_Q;
	model_q->g0_traffic=0.0;
	model_q->g0_valid=0.0;
	model_q->g0_size=0.0;
	model_q->queue_max=10;

	model_q->extra_size=0;
	model_q->extra_traffic=0.;
	model_q->best_extra_size=0;
	model_q->best_extra_traffic=0.;

	mmodel->lba_sampling_ratio=100;
	mmodel->interval_unit_size = ssd_spec->PPS/4;
	mmodel->model_idx = type;
	if (write_size == 0)
		mmodel->time_window=44000*1024/4/mmodel->interval_unit_size*1024; //800GB/(PPB*BPS)
	else
		mmodel->time_window = write_size*1024*1024/4/mmodel->interval_unit_size;

	mmodel->desired_tw = mmodel->time_window;
	mmodel->real_tw = mmodel->time_window*mmodel->interval_unit_size;
	if (type==0) mmodel->model_on=true;
	else mmodel->model_on=false;
	//checking max entry
	mmodel->entry_num = mmodel->real_tw/mmodel->interval_unit_size;
	/*to check first interval
	 * fnumber: # of interval to check
	 * finterval: index of interval to check
	 */
	mmodel->fnumber=10;
	uint32_t finterval[10] = {100, 105, 110, 115, 120, 125, 130, 135, 140, 145};
	for (int i=0;i<mmodel->fnumber;i++) {
		finterval[i] = mmodel->time_window/(mmodel->fnumber+1)*(i+1);
	}
	mmodel->checking_first_interval = (uint32_t*)calloc(mmodel->fnumber, sizeof(uint32_t));
	memcpy(mmodel->checking_first_interval, finterval, sizeof(uint32_t)*mmodel->fnumber);
	mmodel->first_count = 0;
	mmodel->first_interval=0;

	mmodel->one_op=0;
	mmodel->hot_req_count =0;
	mmodel->tot_req_count =0;
	mmodel->hot_lba_count =0;
	mmodel->tot_lba_count =0;
	mmodel->modeling_num = -1;
	mmodel->no_access_lba = 0;
	mmodel->err_cnt=0;
	mmodel->tmp_err_cnt=0;

	/* real or not */
	mmodel->mm_time = (mtime*)calloc(1, sizeof(mtime));
	mmodel->mm_time->is_real=REAL_BENCH;
	
	//TODO have to set interval_unit seg to block
	mmodel->mm_time->interval_unit_size = mmodel->interval_unit_size;	
	mmodel->mm_time->current_time=0;
	mmodel->mm_time->request_time=0;
	mmodel->mm_time->load_timing=0;
	if (mmodel->model_idx == 1) mmodel->mm_time->load_timing = ssd_spec->LBANUM;
	mmodel->mm_time->extra_load = ssd_spec->LBANUM;
	printf("*** MINIATURE MODEL SETTINGS ***\n- lba sampling ratio: %d\n", mmodel->lba_sampling_ratio);
	printf("- interval unit size: %.3f segments (%d pages)\n", 
			(double)mmodel->interval_unit_size/(double)ssd_spec->PPS, mmodel->interval_unit_size);
	printf("- time window: %.2fGB (%llu units)\n", 
			mmodel->time_window/1024./1024.*mmodel->interval_unit_size*4, mmodel->time_window);
	printf("- real mode: %s \n", mmodel->mm_time->is_real?"true":"false");
	printf("- RESIZING: %s \n", RESIZING? "true":"false");
	printf("- extra load timing: %lluGB */\n\n", mmodel->mm_time->extra_load/1024/1024*4);
	
	mmodel->time_stamp = (unsigned long long*)malloc(sizeof(unsigned long long)*ssd_spec->LBANUM/mmodel->lba_sampling_ratio);
	for (int i=0; i<ssd_spec->LBANUM/mmodel->lba_sampling_ratio; i++){
		mmodel->time_stamp[i]=UINT_MAX;
	}
	
	//TODO make max entry (for too many time window)
	mmodel->model_count = (unsigned long long*)calloc(mmodel->entry_num, sizeof(unsigned long long));

	if (type == 0) {
		ginfo = (G_INFO*)malloc(sizeof(G_INFO));
		ginfo->model_ver=-1;
		ginfo->valid=false;
		ginfo->gnum = 0;
		ginfo->commit_g = 19;
		ginfo->vr = (double*)calloc(10, sizeof(double));
		ginfo->commit_vr = (double*)calloc(10, sizeof(double));
		ginfo->gsize = (uint32_t*)calloc(10, sizeof(uint32_t));
		ginfo->WAF = 100.0;
	}
}

/* initialize structure for making miniature model
 */
void model_initialize(mini_model *mmodel, bool on) {
	for (int i=0;i<ssd_spec->LBANUM/mmodel->lba_sampling_ratio; i++) {
		mmodel->time_stamp[i] = UINT_MAX;
	}
	for (int i=0;i<mmodel->entry_num; i++) {
		mmodel->model_count[i] = 0;
	}

	mmodel->first_count=0;
	mmodel->first_interval=0;
	mmodel->mm_time->current_time=0;
	mmodel->mm_time->request_time=0;

	mmodel->one_op = 0;
	mmodel->hot_req_count=0;
	mmodel->tot_req_count=0;
	mmodel->hot_lba_count=0;
	mmodel->tot_lba_count=0;
	mmodel->total_count=0;
	mmodel->no_access_lba = 0;
	mmodel->tmp_err_cnt=0;
	//ginfo->valid=false;
	mmodel->model_on=on;	
}

void model_resizing_timewindow(mini_model *mmodel, unsigned long long timewindow) {
	mmodel->time_window=timewindow;
	mmodel->real_tw = mmodel->time_window*mmodel->interval_unit_size;
	mmodel->entry_num = timewindow;
	free(mmodel->model_count);
	mmodel->model_count = (unsigned long long*)calloc(mmodel->entry_num, sizeof(unsigned long long));
}

/* update time */
void time_managing(mini_model *mmodel, char mode) {
	if (mmodel->mm_time->current_time >= mmodel->time_window) return;
	mmodel->mm_time->request_time++;
	if (mmodel->mm_time->request_time == mmodel->mm_time->interval_unit_size) {
		mmodel->mm_time->current_time++;
		mmodel->mm_time->request_time=0;
	}
}

/* checking time window
 * if return is 0, model is not ready
 * if return is 1, run model
 */
extern bool run_signal;
int check_time_window(mini_model *mmodel, char mode) {
	/* time window end */
	if (mmodel->mm_time->current_time == mmodel->time_window) {
		return 0;
	}

	if (mode == WRITE) {
		/* time update */
		if (mmodel->mm_time->is_real) {
			/* need to check if this is after LOAD_END */
			if (run_signal) {
				if (mmodel->mm_time->load_timing < mmodel->mm_time->extra_load) mmodel->mm_time->load_timing++;
				else time_managing(mmodel, mode);
				
			}
		} else {
			if (mmodel->mm_time->load_timing < (unsigned long)ssd_spec->LBANUM+mmodel->mm_time->extra_load) {
				mmodel->mm_time->load_timing++;
				
			}
			else {
				/* after sequential write */
				time_managing(mmodel,mode);
			}
		}
	}

	/* time window */
	if ((mmodel->mm_time->current_time == 0) && (mmodel->mm_time->request_time == 0)) return 0;
	return 1;
}

/* checking lba's interval */
int check_interval(mini_model *mmodel, uint32_t lba, char mode) {
	/* fixing first interval's count
	 * in this time, we only check one interval unit */
	if (mode == WRITE) {
		
		if (mmodel->model_idx == 0) {
			for (int i=0; i<mmodel->fnumber; i++) {
				if (mmodel->mm_time->current_time < mmodel->checking_first_interval[i]) break;
				if (mmodel->mm_time->current_time == mmodel->checking_first_interval[i]) check_first_interval(mmodel,lba, mode);
			}
		}
		
		
		if (mmodel->mm_time->current_time >= mmodel->time_window) {
			/* collecting count is done
			 * now we have to update last interval count, and make optimal group configuration
			 */
			mmodel->model_on=false;
			mmodel->modeling_num += 1;
			int status = pthread_create(&mmodel->thread_id, NULL, making_group_configuration3, (void*)mmodel);
			if (status != 0) {
				perror("miniature model can't make thread\n");
				abort();
			}
			return 0;
		}
	}
	/* sampling */
	if (lba%mmodel->lba_sampling_ratio) return 0;
	lba = lba/mmodel->lba_sampling_ratio;
	update_count(mmodel,lba, mode);
	return 1;
}

/* updating count by lba */
int update_count(mini_model *mmodel, uint32_t lba, char mode) {
	//TODO
	//check if run correctly
	unsigned long long cur_interval=0;
	unsigned long long tmp_stamp=0;
	unsigned long long tmp_time = mmodel->mm_time->current_time*mmodel->mm_time->interval_unit_size+mmodel->mm_time->request_time;
	
	if (mmodel->time_stamp[lba] == UINT_MAX) {
	       //first accessed
		if (mode == WRITE) {
			mmodel->time_stamp[lba] = tmp_time;
			mmodel->time_stamp[lba] = mmodel->time_stamp[lba] << 1;
			//first access flag
			mmodel->time_stamp[lba]++;
			if (DES_FIX==0) mmodel->tot_req_count ++;
		}
	} else if (mmodel->time_stamp[lba] == UINT_MAX-1) {
		//deleted before this
		if (mode == WRITE) {
			mmodel->time_stamp[lba] = tmp_time;
			mmodel-> time_stamp[lba] = mmodel->time_stamp[lba]<< 1;
			if (DES_FIX==0) mmodel->tot_req_count ++;
		}
	} else if ((mmodel->time_stamp[lba] & 1) == 1) {
		//second accessed
		tmp_stamp = mmodel->time_stamp[lba] >> 1;
		cur_interval = (tmp_time-tmp_stamp)/mmodel->mm_time->interval_unit_size;
		if (cur_interval >= mmodel->entry_num) cur_interval = mmodel->entry_num-1;
		//add 2
		if (cur_interval < 0) {
			abort();
		}
		if (mode == WRITE) {
			mmodel->time_stamp[lba] = tmp_time;
			mmodel->time_stamp[lba] = mmodel->time_stamp[lba]<< 1;
			if (ADD_TWO) mmodel->model_count[cur_interval] += 2;
			else mmodel->model_count[cur_interval] += 1;
		} else {
			mmodel->time_stamp[lba] = UINT_MAX-1;
			mmodel->model_count[cur_interval] ++;
		}
	} else if ((mmodel->time_stamp[lba] & 1) == 0) {
		//or else
		tmp_stamp = mmodel->time_stamp[lba] >> 1;
		cur_interval = (tmp_time-tmp_stamp)/mmodel->mm_time->interval_unit_size;
		if (cur_interval >= mmodel->entry_num) cur_interval = mmodel->entry_num-1;
		if (cur_interval < 0) {
			abort();
		}

		mmodel->model_count[cur_interval]++;
		if (mode == WRITE) {
			mmodel->time_stamp[lba] = tmp_time;
			mmodel->time_stamp[lba] = mmodel->time_stamp[lba]<< 1;
		} else {
			mmodel->time_stamp[lba] = UINT_MAX-1;
		}
	} else {
		abort();
	}
	return 0;
}

double valid_thresh=0.1;
double traffic_thresh=0.15;
bool real_flag=false;
/*complete miniature model and make optimal group configuration using markov chain */
void *making_group_configuration3(void *arg) {
	/* resizing interval count using time stamp (first group, last group, other groups) */
	mini_model *mmodel = (mini_model*)arg;
	mmodel->total_count = 0;
	mmodel->no_access_lba = 0;
	int fixed_G0_size=0;

	//calculate traffic, valid ratio, size of group 0 (HOT)
	model_q->g0_traffic=0;
	model_q->g0_valid=0;
	model_q->g0_size=0;

	model_q->extra_size=0;
	model_q->extra_traffic=0.0;
	model_q->best_extra_size=0;
	model_q->best_extra_traffic=0.0;
	model_q->calc_traffic=0.;
	model_q->calc_size=0;

	for (int i = 0; i < ssd_spec->SEGNUM; i++){
		if (ssd->gnum_info[i] >= ssd_spec->naive_start && ssd->irtable[i] == 0) {
			mmodel->no_access_lba += ssd_spec->PPS;
		}
	}
	mmodel->no_access_lba -= ssd_spec->PPS*ssd->TOTAL_GNUM[ssd_spec->GROUPNUM-1]*0.05;
	if (ssd_spec->mida_on == 1) {
		for (int i = 0; i < 4; i++)
		mmodel->no_access_lba -= ssd_spec->PPS*ssd->TOTAL_GNUM[ssd_spec->naive_start+i]*0.05;
	}

	int qsize = model_q->g0_traffic_queue.size();

	double t=0.0;
	for (int i=0;i<qsize;i++) {
		//model_q->g0_traffic += dqueue_pop_back(&model_q->g0_traffic_queue);
		t = dqueue_pop_back(&model_q->g0_traffic_queue);
		model_q->g0_traffic+=t;
	}
	model_q->g0_traffic = model_q->g0_traffic / (double)qsize;
	if (qsize == 0) {
		model_q->g0_traffic=0.01;
	}

	qsize = model_q->g0_valid_queue.size();
	for (int i=0;i<qsize;i++) model_q->g0_valid += dqueue_pop_back(&model_q->g0_valid_queue);
	model_q->g0_valid = model_q->g0_valid / (double)qsize;
	if (qsize==0) {
		model_q->g0_valid=0.05;
	} 
	
	qsize = model_q->g0_size_queue.size();
	for (int i=0;i<qsize;i++) model_q->g0_size += dqueue_pop_back(&model_q->g0_size_queue);
	model_q->g0_size = model_q->g0_size / (double)qsize;
	model_q->g0_size++; //for active segment
	model_q->g0_size = floor(model_q->g0_size+0.5);
	if (qsize==0) {
		model_q->g0_size=1.0;
	} 
	
	if (model_q->g0_size > 1.0) {
		if (model_q->g0_valid > 0.18) {
			model_resizing_timewindow(mmodel, mmodel->time_window*2);
			model_initialize(mmodel, true);
			return (void*)0;
		}

	} else {
		if (model_q->g0_valid > 0.3) {
			model_resizing_timewindow(mmodel, mmodel->time_window*2);
			model_initialize(mmodel, true);
			return (void*)0;
		}
	}
	if (ssd->TOTAL_GNUM[ssd_spec->GROUPNUM-1]<1) {
		model_resizing_timewindow(mmodel, mmodel->time_window*2);
		model_initialize(mmodel, true);
		return (void*)0;
	}
	
	mmodel->total_count = resizing_model(mmodel);
	
	
	/*making group configuration */
	// there is an opt group config and a valid ratio list for every #ofgroups
	uint32_t **opt_config = (uint32_t**)malloc(sizeof(uint32_t*)*10);
	double **opt_valid_ratio_list = (double**)malloc(10*sizeof(double*));
	for (int i=0;i<10;i++) {
		opt_config[i] = (uint32_t*)calloc(10, sizeof(uint32_t)); //tmp group configuration
		opt_valid_ratio_list[i] = (double*)calloc(10, sizeof(double));
	}
		
	fixed_G0_size = (int)model_q->g0_size;
	//fixed_G0_size = hotfilter->G0_size;
	if (fixed_G0_size == 0) {
		abort();
	}
	opt_config[0][0] = ssd_spec->SEGNUM-ssd_spec->FREENUM;
	opt_valid_ratio_list[0][0] = 1.0;
	double opt_traffic=0.0;
	double opt_waf[10] = {10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0};
	
	//tmp config and tmp valid ratio list
	uint32_t *group_config = (uint32_t*)calloc(10, sizeof(uint32_t));
	double *valid_ratio_list;

	uint32_t gnum=0;
	//save final config and valid ratio list (optimal in every #ofgroups)
	double final_waf = 100.0;
	uint32_t final_gnum=0;
	uint32_t *final_config;
	double *final_valid_ratio_list;
	double final_traffic=0.0;

	//TODO calculate valid ratio of (# of groups 2: HOT, G1)
	mmodel->one_op = 1; // printing option for one group valid ratio predict
	memcpy(group_config, opt_config[0], sizeof(uint32_t)*10);
	
	int g0_desig_size=1;
	double g0_one_traffic=0.;
	double hot_thresh = 0.05;
	if (model_q->g0_valid > valid_thresh || model_q->g0_traffic < traffic_thresh) {
		real_flag = true;
		g0_desig_size=(int)model_q->g0_size;
		group_config[0] -= g0_desig_size;
		model_q->calc_unit=0;
		valid_ratio_list=two_valid_ratio_predictor(mmodel, group_config, 1, mmodel->total_count, g0_desig_size);
		if (valid_ratio_list==NULL) {
		       model_initialize(mmodel, true);
		       return (void*)0;
		}	       
		double predicted_WAF = two_WAF_predictor(mmodel, valid_ratio_list, 1);
		if (predicted_WAF > 100) {
			model_initialize(mmodel, true);
			return (void*)0;
		}
		opt_waf[0] = predicted_WAF;
		opt_config[0][0]=group_config[0];
		opt_valid_ratio_list[0][0]=valid_ratio_list[0];
		model_q->best_extra_size=g0_desig_size;
		model_q->best_extra_traffic=model_q->calc_traffic;
		opt_traffic=model_q->g0_traffic;
	} else {
		real_flag=false;
		for (int i=1;i<30;i++) {
			//model_q->extra_size = i*mmodel->mm_time->interval_unit_size;
			//if (i!=0) group_config[0] -= mmodel->mm_time->interval_unit_size;
			if (i <= model_q->calc_size) continue;
			group_config[0] = ssd_spec->SEGNUM-ssd_spec->FREENUM;
			group_config[0] -= i;
			model_q->calc_unit=0;
			valid_ratio_list = two_valid_ratio_predictor(mmodel,group_config, 1, mmodel->total_count, i);
			if (valid_ratio_list==NULL) break;
			double predicted_WAF = two_WAF_predictor(mmodel,valid_ratio_list, 1);
			//for the case that extra_cnt is not same as real size
			if (i==1) {
				g0_one_traffic=model_q->calc_traffic;
				//TODO decide when hot group will be fixed
				if (HOT_POSITION == 1) {
					opt_waf[0]=predicted_WAF;
					opt_config[0][0]=group_config[0];
					opt_valid_ratio_list[0][0]=valid_ratio_list[0];
					model_q->best_extra_size=model_q->calc_size;
					model_q->best_extra_traffic=g0_one_traffic;
					model_q->best_extra_unit=model_q->calc_unit;
					free(valid_ratio_list);
					break;
				}
			}
			if (predicted_WAF < opt_waf[0]-hot_thresh) {
				model_q->best_extra_size=model_q->calc_size;
				model_q->best_extra_traffic=model_q->calc_traffic;
				model_q->best_extra_unit=model_q->calc_unit;
				opt_waf[0]=predicted_WAF;
				opt_config[0][0] = group_config[0];
				opt_valid_ratio_list[0][0] = valid_ratio_list[0];
			} else if (i==0) {
				free(valid_ratio_list);
				continue;
			} else {
				free(valid_ratio_list);
				break;
			}
		}
		if (HOT_POSITION == 0) {
			model_q->calc_traffic = model_q->best_extra_traffic;
			g0_desig_size=model_q->best_extra_size;
			model_q->calc_unit = model_q->best_extra_unit;
			opt_traffic = model_q->calc_traffic;
		} else {
			model_q->calc_traffic = g0_one_traffic;
			g0_desig_size=1;
			model_q->calc_unit = model_q->best_extra_unit;
			opt_traffic = model_q->calc_traffic;
		}
	}


	int group_add_flag=0;
	uint32_t max_groupsize=opt_config[0][0]-50;
	uint32_t tmp_maxsize = max_groupsize;
	/* making group configurations */
	for (int i=0; i<8; i++) {
		if ((mmodel->model_idx == 0) && (i==8)) break;
		gnum = i+2;
		group_add_flag=0;
		memcpy(group_config, opt_config[i], sizeof(uint32_t)*10);

		group_config[i+1] = group_config[i];
		group_config[i] = 1;
		max_groupsize = tmp_maxsize;
		for (int j=1;j<max_groupsize; j++) {
			group_config[i]=j;
			group_config[i+1] -= 1;
			if (group_config[i+1]==0) {
				abort();
				break;
			}
			/* predict valid ratio for the group configuration */
			valid_ratio_list = two_valid_ratio_predictor(mmodel,group_config, gnum, mmodel->total_count, g0_desig_size);
			if (valid_ratio_list == NULL) {
			} else {
			/* calculate WAF using markov chain */
				double predicted_WAF = two_WAF_predictor(mmodel,valid_ratio_list, gnum);
				if (opt_waf[i+1] > predicted_WAF) {
					//opt_gnum = gnum;
					opt_waf[i+1] = predicted_WAF;
					memcpy(opt_config[i+1], group_config, sizeof(uint32_t)*10);
					memcpy(opt_valid_ratio_list[i+1], valid_ratio_list, sizeof(double)*10);
					group_add_flag=1;
					tmp_maxsize = max_groupsize - opt_config[i+1][i];
				}
				free(valid_ratio_list);
			}
		}
		memcpy(group_config, opt_config[i+1], sizeof(uint32_t)*10);
	}

	for (int i=0;i<10;i++) {
		if (final_waf > opt_waf[i]) {
			final_gnum = i+1;
			final_config = opt_config[i];
			final_valid_ratio_list = opt_valid_ratio_list[i];
			final_waf = opt_waf[i];
		}
		else {
			free(opt_config[i]);
			free(opt_valid_ratio_list[i]);
		}
	}

	int initial_last_size=0;
	model_q->calc_size=0;
	/* group configuration phase 2*/
	for (int n=0;n<4;n++) {

		if (final_waf == 100) break;
		if (real_flag==false) {
			if (HOT_POSITION==0) {
				memcpy(group_config, final_config, sizeof(uint32_t)*10);
				gnum=final_gnum;
				
				group_config[gnum-1] += g0_desig_size;
				initial_last_size = group_config[gnum-1];
				model_q->calc_size=0;	
	
				for (int i=1;i<30;i++) {
					if (i <= model_q->calc_size) continue;
					group_config[gnum-1] = initial_last_size;
					group_config[gnum-1] -= i;
					model_q->calc_unit=0;
					valid_ratio_list=two_valid_ratio_predictor(mmodel, group_config, gnum, mmodel->total_count, i);
					if (valid_ratio_list == NULL) {
						break;
					}
					double predicted_WAF = two_WAF_predictor(mmodel, valid_ratio_list, gnum);
	
					if (predicted_WAF < final_waf) {
						model_q->best_extra_size=model_q->calc_size;
						model_q->best_extra_traffic=model_q->calc_traffic;
						model_q->best_extra_unit = model_q->calc_unit;
						final_waf=predicted_WAF;
						memcpy(final_config, group_config, sizeof(uint32_t)*10);
						memcpy(final_valid_ratio_list, valid_ratio_list, sizeof(double)*10);
						g0_desig_size=i;
						opt_traffic = model_q->calc_traffic;
					} else if (i==1) continue;
					else continue;
					free(valid_ratio_list);
				}
			}
		}

		gnum = final_gnum;
		memcpy(group_config, final_config, sizeof(uint32_t)*10);
		model_q->calc_unit = model_q->best_extra_unit;
		for (int i=0;i<8;i++) {
			if ((mmodel->model_idx == 0) && (i==9)) break;
			if (i < gnum-1) {
				group_config[gnum-1] += group_config[i];
			} else {
				group_config[i+1] = group_config[i];
				gnum++;
			}
			max_groupsize = group_config[gnum-1]-50;
			group_add_flag=0;
			group_config[i] = 0;
			for (int j=1;j<max_groupsize; j++) {
				group_config[i] = j;
				group_config[gnum-1] -= 1;
				valid_ratio_list = two_valid_ratio_predictor(mmodel,group_config, gnum, mmodel->total_count, g0_desig_size);
				if (valid_ratio_list == NULL) {
					break;
				}
				double predicted_WAF = two_WAF_predictor(mmodel,valid_ratio_list, gnum);
				if (final_waf > predicted_WAF) {
					final_gnum = gnum;
					final_waf = predicted_WAF;
					memcpy(final_config, group_config, sizeof(uint32_t)*10);
					memcpy(final_valid_ratio_list, valid_ratio_list, sizeof(double)*10);
					group_add_flag=1;
				}
				free(valid_ratio_list);
			}
			memcpy(group_config, final_config, sizeof(uint32_t)*10);
			if (group_add_flag==0) {
				if (i >= final_gnum-1) break;
			}
			
		}
	}
	
	if (HOT_POSITION && (model_q->g0_valid <= valid_thresh) && (model_q->g0_traffic >= traffic_thresh)) {
	
		memcpy(group_config, final_config, sizeof(uint32_t)*10);
		gnum = final_gnum;

		for (int i=2;i<30;i++) {
			group_config[gnum-1]--;
			valid_ratio_list =two_valid_ratio_predictor(mmodel,group_config, gnum, mmodel->total_count, i);
			if (valid_ratio_list==NULL) {
				break;
			}
			double predicted_WAF = two_WAF_predictor(mmodel,valid_ratio_list, gnum);
			if (predicted_WAF < final_waf) {
				model_q->best_extra_size=i;
				model_q->best_extra_traffic=model_q->calc_traffic;
				final_waf=predicted_WAF;
				memcpy(final_config, group_config, sizeof(uint32_t)*10);
				memcpy(final_valid_ratio_list, valid_ratio_list, sizeof(double)*10);
				opt_traffic = model_q->calc_traffic;
			} else {
				continue;
			}	
		}
	}


	if (mmodel->model_idx==0) {
		final_gnum++;
		for (int i=final_gnum-2;i>=0;i--) {
			final_config[i+1] = final_config[i];
			final_valid_ratio_list[i+1] = final_valid_ratio_list[i];
		}
		final_config[0] = model_q->best_extra_size;
		final_valid_ratio_list[0]=0.0;
		final_traffic = opt_traffic;
	}

	//set group information
	if (final_gnum > 2) {
		ginfo->model_ver = mmodel->model_idx;
		memcpy(ginfo->gsize, final_config, sizeof(uint32_t)*10);
		memcpy(ginfo->vr, final_valid_ratio_list, sizeof(double)*10);
		ginfo->gnum = final_gnum;
		ginfo->WAF = final_waf;
		ginfo->g0_traffic = final_traffic;
		print_config_new(ginfo->gnum, ginfo->gsize, ginfo->WAF, ginfo->vr, ginfo->g0_traffic);
		ginfo->valid=true;
	} 

	return (void*)0;
}

/* print final results by miniature model */
void print_config_new(int gnum, uint32_t *opt_config, double opt_waf, double *opt_valid_ratio_list, double tr) {
	printf("\n*****MODEL PREDICTION RESULTS*****\n");
	printf("group number: %d\n", gnum);
	for (int i=0;i<gnum; i++) printf("*group %d: size %d, valid ratio %f\n", i, opt_config[i], opt_valid_ratio_list[i]);
	printf("calculated WAF: %f\n", opt_waf);
	if (real_flag == true) printf("Used traffic (REAL) : %.3f\n", tr);
	else printf("Used traffic (CALC) : %.3f\n", tr);
	printf("************************************\n\n");
}

void print_config_into_log(int gnum, uint32_t *opt_config, double opt_waf, double *opt_valid_ratio_list) {
	char l[128] = "logging/log_";
	strcat(l, workload_name);
	FILE* mFile = fopen(l, "a");
	fprintf(mFile, "*****MODEL PREDICTION RESULTS*****\n");
	fprintf(mFile, "group number: %d\n", gnum);
	for (int i=0;i<gnum; i++) fprintf(mFile, "*group %d: size %d, valid ratio %f\n", i, opt_config[i], opt_valid_ratio_list[i]);
	fprintf(mFile, "calculated WAF: %f\n", opt_waf);
	fprintf(mFile, "************************************\n");
	fclose(mFile);
}


/* predict WAF by group's valid ratios */
double two_WAF_predictor(mini_model *mmodel,double *valid_ratio_list, int group_num) {
	double **valid_matrix = (double**)malloc(sizeof(double*)*(group_num+2));
	double *valid_ratio_list2 = (double*)malloc(sizeof(double)*(group_num+2));
	valid_ratio_list2[0] = 1.0;
	//TODO check g0 valid ratio
	if (real_flag == true) valid_ratio_list2[1] = model_q->g0_valid;
	else valid_ratio_list2[1] = 0.0;
	memcpy(&valid_ratio_list2[2], valid_ratio_list, sizeof(double)*group_num);

	for (int i=0;i<(group_num+2); i++) {
		double *tmp = (double*)calloc(group_num+2, sizeof(double));
		tmp[0] = 1.0-valid_ratio_list2[i];
		if (i==0) {
			if (real_flag==true) tmp[1] = model_q->g0_traffic;
			else tmp[1] = model_q->calc_traffic;
			tmp[2] = 1.0 - tmp[1];
		}
		else if (i==group_num+1) tmp[i] = valid_ratio_list2[i];
		else tmp[i+1] = valid_ratio_list2[i];
		valid_matrix[i]=tmp;
	}	
	
	double *base = (double*)calloc(group_num+2, sizeof(double));
	base[0]=ssd_spec->PGNUM*0.1;
	base[1]=ssd_spec->PGNUM*0.9;
	double *result = (double*)calloc(group_num+2, sizeof(double));
	memcpy(result, base, sizeof(double)*(group_num+2));
	uint32_t past=0;
	uint32_t write=0;
	double total_wr=0.0;
	double tmp_tot_wr=0.0;

	for (int i=0;i<2000;i++) {
		tmp_tot_wr=0;
		memcpy(base, result, sizeof(double)*(group_num+2));
		memset(result, 0, sizeof(double)*(group_num+2));
		//dot operation
		for (int j=0;j<(group_num+2);j++) {
			for (int k=0;k<(group_num+2);k++) {
				result[j] += base[k]*valid_matrix[k][j];
			}
		}

		if (group_num == 1) {
			total_wr += result[2]-result[0]*valid_ratio_list2[1];
			tmp_tot_wr += result[2]-result[0]*valid_ratio_list2[1];
		} else {
			for (int j=0;j<(group_num-1);j++) {
				total_wr += result[j+3];
				tmp_tot_wr += result[j+3];
			}
			total_wr += result[1]*valid_ratio_list2[1];
			tmp_tot_wr += result[1]*valid_ratio_list2[1];
		}
		past = result[0];
		write += past;
	}
	double WAF = (tmp_tot_wr+(double)past)/(double)past;

	/* free memory */
	free(result);
	free(base);
	for (int i=0;i<(group_num+2);i++) free(valid_matrix[i]);
	free(valid_ratio_list2);
	free(valid_matrix);


	return WAF;
}


/* resize model's interval counts (first group, last group, other groups) */
unsigned long long resizing_model(mini_model *mmodel) {

	/* first interval unit */
	char m1[128] = "modeling/modeling_";
	char m2[128] = "modeling/models_";
	//char m2[128] = "modeling/modeling_prev_";
	strcat(m1, workload_name);
	
	strcat(m2, workload_name);
	strcat(m2, "/");
	strcat(m2, "model_");
	char cmn[128];
	sprintf(cmn, "%d", mmodel->modeling_num);
	strcat(m2, cmn);

	strcat(m2, workload_name);

	/* last interval unit */
	//mmodel->model_count[mmodel->time_window-1] = utilization/mmodel->lba_sampling_ratio;
	for (int i=0;i<ssd_spec->LBANUM/mmodel->lba_sampling_ratio; i++) {
		if ((mmodel->time_stamp[i] != UINT_MAX) && ((mmodel->time_stamp[i]&1) == 1) && (mmodel->time_stamp[i] != UINT_MAX-1)) {
			mmodel->model_count[mmodel->entry_num-1]++;
		}
	}
	unsigned long long tot=0;
	/* first interval unit */
	if (FIRST_NEW == 0) {
		if ((mmodel->lba_sampling_ratio != 1) && (mmodel->fnumber != 1)) {
			mmodel->model_count[0] = mmodel->first_count * mmodel->time_window /mmodel->fnumber/ mmodel->lba_sampling_ratio;
		}
	} else {
		mmodel->model_count[0] = mmodel->first_interval/mmodel->lba_sampling_ratio;
		unsigned long long tmp = mmodel->first_count * mmodel->time_window /mmodel->fnumber/ mmodel->lba_sampling_ratio;
	}
	tot += mmodel->model_count[0];
	
	double left_time = 0.0;
	/* other interval units */
	for (int i=1;i<mmodel->entry_num;i++) {
		if (RESIZING == 1) {
			left_time = ((double)(mmodel->time_window%i*2+i)/2.0)/(double)i;
			mmodel->model_count[i] = mmodel->model_count[i] + (int)((left_time*mmodel->model_count[i])/((double)mmodel->time_window/(double)(i)));
		}
		tot += mmodel->model_count[i];
	}

	return tot;
}

double FIFO_predict(double op){
	return op/(op+lambert_w0(-op*exp(-op)));	
}
double None_FIFO_predict(mini_model *mmodel, double op){
	if(op<=1){
		return 0.;
	}
	double r = (double) mmodel->hot_req_count / (double) mmodel->tot_req_count;
	double f = (double) mmodel->hot_lba_count / (double) mmodel->tot_lba_count;
	
	if (mmodel->one_op == 1){
		mmodel->one_op = 0;
	}
	double inc_waf = 1.;
	double test_val = 1.;
	int inf_flag = 0;
	while(test_val > 0){
		test_val = 1 + r/double(exp(double(r*op)/double(f*inc_waf))-1) + (1-r)/double(exp(double((1-r)*op)/double(double(1-f)*inc_waf))-1)-inc_waf;
		inc_waf += 1;
	}
	inc_waf -= 2;
	test_val = 1.;
	while(test_val > 0){
		test_val = 1 + r/double(exp(double(r*op)/double(f*inc_waf))-1) + (1-r)/double(exp(double((1-r)*op)/double(double(1-f)*inc_waf))-1)-inc_waf;
		inc_waf += 0.1;
	}
	inc_waf -= 0.2;
	test_val = 1.;
	while(test_val > 0){
		test_val = 1 + r/double(exp(double(r*op)/double(f*inc_waf))-1) + (1-r)/double(exp(double((1-r)*op)/double(double(1-f)*inc_waf))-1)-inc_waf;
		inc_waf += 0.01;
	}
	inc_waf -= 0.02;
	test_val = 1.;
	while(test_val > 0){
		test_val = 1 + r/double(exp(double(r*op)/double(f*inc_waf))-1) + (1-r)/double(exp(double((1-r)*op)/double(double(1-f)*inc_waf))-1)-inc_waf;
		inc_waf += 0.001;
	}
	inc_waf -= 0.001;
	
	return inc_waf;
}


double one_group_predictor(mini_model *mmodel, uint32_t *group_config, uint32_t group_num, unsigned long long tot_cnt) {	//input values are same as vlid_ratio_predictor, return valid ratio
	double last_vr;
	double last_group_vp = (double)utilization;  //global var
	double last_waf = None_FIFO_predict(mmodel,(double)group_config[group_num-1]*ssd_spec->PPS/(double)last_group_vp);	//None_FIFO_predict()'s input is op 
	if (DES == 1){	//FIFO
		last_vr = 1.-1./last_waf;
		return last_waf;
	}else if (DES == 2){	//Greedy
		last_waf = None_FIFO_predict(mmodel,(1+1/double(2*ssd_spec->PPS))*(double)group_config[group_num-1]*ssd_spec->PPS/(double)last_group_vp)/(1+1/double(2*ssd_spec->PPS));	//Changed op is input
		last_vr = 1.-1./last_waf;
		return last_waf;
	}
	else{
		return -1;
	}
	
}

double *two_valid_ratio_predictor(mini_model *mmodel, uint32_t *group_config, uint32_t group_num, unsigned long long tot_cnt, int extra_cnt) {
	double *valid_ratio_list = (double*)calloc(10, sizeof(double));
	unsigned long long *invalid_cnt_list = (unsigned long long*)calloc(10, sizeof(unsigned long long));
	
	uint32_t groups_seg_num=0;
	for (int i=0;i<group_num-1;i++) groups_seg_num += group_config[i];
	double mig_rate = 1.0; //for calcuting group_time
	double ref_rate = 1.0;
	double tmp_vr=0.0;
	double tmp_tr=0.0;
	double valid_ratio=0.0;

	unsigned long long invalid_cnt = 0;
	int extra_g0_cnt=0;
	double valid_cnt = 0.0;
	double global_valid_cnt = 0.0;

	double cur_interval=0.0;
	double group_time = 0.0;
	double n_group_time=0.0; //ROUND(group_time)
	//Jeeyun check segment time
	double one_seg_time = 0.0;
	double one_invalid_cnt = 0.0;

	double past_group_time = 0.0;

	unsigned long long lindex=0;

	double evg_lifespan = 0.0;
	unsigned long long lg_req = 0;
	unsigned long long lifespan_num = 0;

	double valid_g0=0.;
	double invalid_g0=0.;
	int real_g0_entry=0;
	double calculated_traffic=0.;

	model_q->extra_traffic=0.;
	model_q->extra_size=0;

	//for the case that the group 0 size is too small
	if ((double)group_config[0] < (double)mmodel->mm_time->interval_unit_size/(double)ssd_spec->PPS) {
		for (int i=0;i<group_num;i++) valid_ratio_list[i] = 1.0;
		free(invalid_cnt_list);
		return valid_ratio_list;
	}

	double extr_size=0.;
	double size_gap=0.;
	double valid_thresh=0.2;
	//calculate the g0 traffic of the designated g0 size
	if (real_flag==false) {
		//can calculate the HOT group size
		if (model_q->calc_unit == 0) {
			//need to calculate the traffic and the unit num
			for (int i=0;i<mmodel->entry_num;i++) {
				invalid_cnt += mmodel->model_count[i];
				tmp_vr = 1.0-(double)invalid_cnt/(double)tot_cnt;
				tmp_tr = (double)invalid_cnt/(double)tot_cnt;
				cur_interval = (double)(i+1)*(double)mmodel->mm_time->interval_unit_size/(double)ssd_spec->PPS;
				extr_size = floor((double)cur_interval*tmp_tr+0.5);
				if (extr_size >= extra_cnt) {
					past_group_time += cur_interval;
					model_q->calc_traffic=tmp_tr;
					tot_cnt -= invalid_cnt;
					if (extr_size > extra_cnt) {
						size_gap = extr_size - extra_cnt;
						group_config[group_num-1] -= size_gap;
					}
					model_q->calc_size=extr_size;
					model_q->calc_unit=i+1;
					break;
				}
			}
		} else {
			//there is traffic and the unit num info
			for (int i=0;i<model_q->calc_unit;i++) {
				invalid_cnt += mmodel->model_count[i];
				tmp_vr = 1.0-(double)invalid_cnt/(double)tot_cnt;
				global_valid_cnt += tmp_vr*mmodel->mm_time->interval_unit_size*ref_rate;
			}
			tmp_tr = (double)invalid_cnt/(double)tot_cnt;
			cur_interval = (double)model_q->calc_unit*
				(double)mmodel->mm_time->interval_unit_size/(double)ssd_spec->PPS;
			extr_size = floor((double)cur_interval*tmp_tr+0.5);
			past_group_time += cur_interval;
			model_q->calc_traffic=tmp_tr;
			tot_cnt -= invalid_cnt;
			model_q->calc_size = extr_size;
		}
	} else {
		//g0 size is fixed
		mig_rate=1.0/(model_q->g0_traffic);
		ref_rate=model_q->g0_traffic;
		group_time = (double)extra_cnt*mig_rate-mig_rate/2.0;
		for (int i=0;i<mmodel->entry_num;i++) {
			cur_interval=(double)i*(double)mmodel->mm_time->interval_unit_size/ssd_spec->PPS;
			if (cur_interval<group_time) {
				invalid_cnt += mmodel->model_count[i];
				tmp_vr = 1.0-(double)(invalid_cnt)/(double)(tot_cnt);
				global_valid_cnt += tmp_vr*ref_rate*mmodel->mm_time->interval_unit_size;
			} else {
				real_g0_entry=i;
				model_q->calc_unit = i;
				break;
			}
		}
		past_group_time += group_time;
		model_q->calc_traffic=(double)invalid_cnt/(double)tot_cnt;
		tot_cnt -= invalid_cnt;
	}


	
	double g1_traffic = 1.0-model_q->calc_traffic;
	if (real_flag==true) g1_traffic = 1.0-model_q->g0_traffic
		+(model_q->g0_traffic*model_q->g0_valid);
	for (int i=0;i<group_num; i++) {
		invalid_cnt = 0;
		valid_cnt = 0.0;
		one_invalid_cnt=0.0;
		/* calculate group_time value 
		 * (increased when group number increase)
		 */
		
		mig_rate = 1.0/g1_traffic;
		ref_rate = g1_traffic;
		for (int j=0;j<i;j++) {
			if (j==group_num-1) {
				double tmp = valid_ratio_list[group_num-2]+valid_ratio_list[group_num-1];
				mig_rate *= (1.0/(double)tmp);
				break;
			} else {
				mig_rate *= (1.0/valid_ratio_list[j]);
				ref_rate *= valid_ratio_list[j];
			}
		}
		group_time = (double)group_config[i]*mig_rate-mig_rate/2;
		//error check: if number of segments of groups over the total time window = time window < device size
		if ((i!=group_num-1) && (((uint32_t)group_time+past_group_time) >= ((mmodel->time_window*mmodel->mm_time->interval_unit_size/ssd_spec->PPS-50)))) {
			if (model_q->g0_valid > valid_thresh && group_num==1) {
				abort();
			}

			free(invalid_cnt_list);
			return NULL;
		}

		/* add invalid counts
		 */
		unsigned long long tmp2=0;
		for (int k=0;k<i;k++) tmp2 += invalid_cnt_list[k];
		for (int j=0;j<mmodel->entry_num;j++) {
			cur_interval = (double)(j)*(double)mmodel->mm_time->interval_unit_size/(double)ssd_spec->PPS;
			if ((cur_interval >= past_group_time) && (cur_interval < group_time+past_group_time)) {
				
				invalid_cnt += mmodel->model_count[j];
				
				if (i < group_num-1) {
					//tmp_vr: expected valid ratio of that group
					tmp_vr = 1.0-(double)(invalid_cnt)/(double)(tot_cnt-tmp2);
					global_valid_cnt += tmp_vr*ref_rate*mmodel->mm_time->interval_unit_size;
				}else if (i==group_num-1) {
				}	
			}
			if ((i == group_num-1) && (cur_interval >= past_group_time)) {
				if ((DES_FIX == 0) && (cur_interval >= group_time+past_group_time)) break;
				lg_req += mmodel->model_count[j]*cur_interval;
				lifespan_num += mmodel->model_count[j];
			}
		}

		//TODO testing extra g0 cnt
		if (i==0) invalid_cnt += extra_g0_cnt;
		// for Desnoyer for last group
		if (i==group_num-1) {
			evg_lifespan = (double) lg_req / (double) lifespan_num ;  //caclulate the average of lifespan in Last group
			for (int j=0;j<mmodel->entry_num;j++) {
				cur_interval = (double)(j)*(double)mmodel->mm_time->interval_unit_size/(double)ssd_spec->PPS;
				if ((cur_interval >= past_group_time)) {
					if ((DES_FIX == 0) && (cur_interval >= group_time+past_group_time)) break;
					mmodel->tot_req_count += mmodel->model_count[j]; // caculate total request count in LAST GROUP
					double lba_num = (double) mmodel->model_count[j] / ((double) ((mmodel->time_window*mmodel->mm_time->interval_unit_size-50)) / (double) cur_interval); //count lba that has lifespan as cur_interval
					mmodel->tot_lba_count += lba_num; // caculate total lba count in LAST GROUP
					if (DES_FIX==1) {
						if(cur_interval <= evg_lifespan*0.7){	//check is this request hot
							mmodel->hot_lba_count += lba_num;
							mmodel->hot_req_count += mmodel->model_count[j];
						}
					} else {
						if ((cur_interval < group_time+past_group_time) && (j != mmodel->entry_num-1)) {
						mmodel->hot_lba_count += lba_num;
						mmodel->hot_req_count += mmodel->model_count[j];
					}
					}
				}
			}
		}
		//update invalid count
		//CHECK invalid cnt is mean value of last one segment's invalid cnt
		invalid_cnt_list[i] = invalid_cnt;
		//calculate valid ratio
		unsigned long long tmp3 = 0;
		for (int j=0;j<i;j++) tmp3 += invalid_cnt_list[j];
		valid_ratio = 1.0-(double)invalid_cnt/(double)(tot_cnt-tmp3);
		//TODO last group is useless, need to change code
		if (i < group_num-1) {
			valid_ratio_list[i] = valid_ratio;
			past_group_time += (group_time);
			//if (i==0) past_group_time++;
		}
	}

	double new_util;
	if ((mmodel->no_access_lba < 0) || (mmodel->no_access_lba < ssd_spec->SEGNUM*ssd_spec->PPS*0.2)){
		mmodel->no_access_lba = 0;
	}

	new_util = (double)utilization - (double)mmodel->no_access_lba;
	double new_group_size = (double)group_config[group_num-1]*ssd_spec->PPS - (double)mmodel->no_access_lba;
	if (new_util < global_valid_cnt) {
		return NULL;
	}
	double last_group_vp = (double)new_util - global_valid_cnt;
	//double last_group_vp = (double)utilization - global_valid_cnt;
	if (last_group_vp < 0) {
		free(invalid_cnt_list);
		return NULL;
	}
	double last_waf = None_FIFO_predict(mmodel,new_group_size/last_group_vp);
	double last_vr;
	last_vr = 1.-1./last_waf;
	
	if (last_vr > 1.0) {
		free(invalid_cnt_list);
		return NULL;
	}
	
	if (DES == 1){ // FIFO
		valid_ratio_list[group_num-1] = last_vr;
		goto END;
	}
	if (DES == 2){ // Greedy
		last_waf = None_FIFO_predict(mmodel,(1+1/double(2*ssd_spec->PPS))*(double)group_config[group_num-1]*ssd_spec->PPS/(double)last_group_vp)/(1+1/double(2*ssd_spec->PPS));
		last_vr = 1.-1./last_waf;
		valid_ratio_list[group_num-1] = last_vr;
		goto END;
	}
		
END:
//	if (valid_ratio_list[group_num-1]<=0) {
	for(int k = 0; k < group_num; k++){
		if (valid_ratio_list[k] <= 0.){	
			free(invalid_cnt_list);
			group_config[group_num-1] += model_q->extra_size;
			return NULL;
		}
	}

	/* free memory */
	free(invalid_cnt_list);
	group_config[group_num-1] += model_q->extra_size;
	return valid_ratio_list;
}

/*initialize data structures for first interval search*/
void initialize_first_interval(mini_model *mmodel) {
	mmodel->updated_lbas = (uint32_t*)malloc(sizeof(uint32_t)*(mmodel->mm_time->interval_unit_size+1));
	for (int i=0;i<mmodel->mm_time->interval_unit_size+1; i++) mmodel->updated_lbas[i]=UINT_MAX;
}

/* fixing first interval's count
 * there is no sampling */
int check_first_interval(mini_model *mmodel, uint32_t lba, char mode) {
	//initialize
	if (mmodel->mm_time->request_time == 0) {
		initialize_first_interval(mmodel);
	}
	//lba search
	int j=0;
	int indx=-1;
	while (mmodel->updated_lbas[j] != UINT_MAX) {
		if (mmodel->updated_lbas[j] == lba) {
			indx = j;
			break;
		}
		j++;
	}
	//lba check
	if (indx != -1) mmodel->first_count++;
	else mmodel->updated_lbas[j] = lba;

	if (mmodel->mm_time->request_time == mmodel->mm_time->interval_unit_size-1) {
		//remove
		remove_first_interval(mmodel);
	}	
	return 1;
}

void remove_first_interval(mini_model *mmodel) {
	int j=0;
	while (mmodel->updated_lbas[j] != UINT_MAX) j++;
	free(mmodel->updated_lbas);
}


void model_destroy(mini_model *mmodel) {
	//free(mmodel->time_stamp);
	free(mmodel->model_count);
	free(mmodel);
}
