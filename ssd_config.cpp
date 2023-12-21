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

#include "ssd_config.h"
#include "queue.h"

extern SSD_SPEC *ssd_spec;
extern struct SSD *ssd;
extern deque<int> fbqueue;

void stats_init(struct STATS *stats){

	stats->cur_req = 0; //current # of write request
	stats->cur_time=0;
	stats->cur_wp = 1; //current # of write page
	stats->copy = 0; //current # of copy page
	stats->tmp_wp = 1; //current # of write page in time window
	stats->tmp_copy = 0; //current # of copy page in time window
	stats->erase = 0; //current # of block erase
	stats->tot_req = 0; //total # of write iequest
	stats->progress = 0; // progress of simulation
	stats->util=0;
	stats->vp=0;
	stats->fourty_waf=0.;
	stats->fourty_wp=1.;
	stats->unknown_req = 0;
	stats->unknown_cnt = 0;
	stats->m_wp = 1;
	stats->m_copy = 0;
	stats->checktime = 0.;
	stats->err_write = 0;
	stats->err_flag = 0;
	stats->commit_bit = 0;
}

int policy_to_int(char *policy){
	if(strcmp(policy, "fifo") == 0) return 0;
	else if(strcmp(policy, "greedy") == 0) return 1;
	else if(strcmp(policy, "cost-benefit") == 0) return 2;
	else return 1;
}

void ssd_init(struct SSD* ssd, class GROUP* group[], int gnum, int vs_policy, int naive_start, int dev_gb, int seg_mb){
	//Spec initialization
	ssd_spec = (struct SSD_SPEC*)malloc(sizeof(struct SSD_SPEC));
	ssd_spec->MAXGNUM = 20; //for group configuration change (not used yet)
	ssd_spec->GROUPNUM = gnum;
	ssd_spec->vs_policy = vs_policy;
	ssd->config_version=0;

	ssd_spec->LOGSIZE = (long)dev_gb*1000*1000*1000; //logical size
	ssd_spec->DEVSIZE = (long)dev_gb*1024*1024*1024; //physical size

	ssd_spec->dev_gb = dev_gb;
	ssd_spec->seg_mb = seg_mb;

	ssd_spec->naive_start = naive_start;
	ssd_spec->mida_on = 1;
	ssd_spec->PGSIZE = 4096; //4KB
	ssd_spec->PPB = 512; // page per block
	ssd_spec->BPS = seg_mb/2; // block per segment to 128 
	ssd_spec->PPS = ssd_spec->PPB*ssd_spec->BPS;
	ssd_spec->FREENUM = 4; //GC trigger point
	ssd_spec->NAIVE_N = 4; 
	ssd_spec->BLKSIZE = ssd_spec->PGSIZE*ssd_spec->PPB;
	ssd_spec->SEGSIZE = ssd_spec->BLKSIZE*ssd_spec->BPS;
	ssd_spec->SEGNUM = (int)(ssd_spec->DEVSIZE/ssd_spec->SEGSIZE);
	ssd_spec->BLKNUM = ssd_spec->SEGNUM*ssd_spec->BPS;
	ssd_spec->PGNUM = ssd_spec->BLKNUM*ssd_spec->PPB;
	ssd_spec->LBANUM = (int)(ssd_spec->LOGSIZE/ssd_spec->PGSIZE);
	//simulation data structure
	ssd->itable = (bool*)malloc(sizeof(bool)*(ssd_spec->PGNUM+1));
	ssd->mtable = (int*)malloc(sizeof(int)*ssd_spec->LBANUM);
	ssd->oob = (struct OOB*)malloc(sizeof(struct OOB)*ssd_spec->PGNUM);

	//malloc data structure
	ssd->irtable = (int*)malloc(sizeof(int)*ssd_spec->SEGNUM);
	ssd->gnum_info = (unsigned char*)malloc(sizeof(unsigned char)*ssd_spec->SEGNUM);
	ssd->fill_info = (bool*)malloc(sizeof(bool)*ssd_spec->SEGNUM);
	ssd->seg_stamp = (double*)malloc(sizeof(double)*ssd_spec->SEGNUM);
	ssd->VPG = (double*)malloc(sizeof(double)*ssd_spec->MAXGNUM);
	ssd->EPG = (double*)malloc(sizeof(double)*ssd_spec->MAXGNUM);
	ssd->DEBUG_VPG = (double*)malloc(sizeof(double)*ssd_spec->MAXGNUM);
	ssd->DEBUG_EPG = (double*)malloc(sizeof(double)*ssd_spec->MAXGNUM);
	ssd->TMP_VPG = (double*)malloc(sizeof(double)*ssd_spec->MAXGNUM);
	ssd->TMP_EPG = (double*)malloc(sizeof(double)*ssd_spec->MAXGNUM);
	ssd->CUR_VPG = (double*)malloc(sizeof(double)*ssd_spec->MAXGNUM);
	ssd->TOTAL_GNUM = (double*)malloc(sizeof(double)*ssd_spec->MAXGNUM);
	ssd->active = (int**) malloc(sizeof(int*)*ssd_spec->MAXGNUM);

	ssd->err_VPG = (double*)malloc(sizeof(double)*ssd_spec->MAXGNUM);
	ssd->err_EPG = (double*)malloc(sizeof(double)*ssd_spec->MAXGNUM);


	for (int i = 0; i < ssd_spec->MAXGNUM; i++){
		ssd->active[i] = (int*)malloc(sizeof(int)*2);
		//first element: idx of active block, second element: appending point on active block
		ssd->active[i][0] = -1;
		ssd->active[i][1] = 0;
	}
	printf("here\n");
	//free queue management
	for(int i = 0; i < ssd_spec->SEGNUM; i++){
		fbqueue.push_front(i);
	}
	printf("here\n");

	//data structure initialization
	memset(ssd->itable, 0, sizeof(bool)*ssd_spec->PGNUM);
	memset(ssd->mtable, -1, sizeof(int)*ssd_spec->LBANUM);
	memset(ssd->irtable, 0, sizeof(int)*ssd_spec->SEGNUM);
	memset(ssd->gnum_info, ssd_spec->MAXGNUM-1, sizeof(unsigned char)*ssd_spec->SEGNUM);
	memset(ssd->fill_info, 0, sizeof(bool)*ssd_spec->SEGNUM);
	memset(ssd->seg_stamp, 0., sizeof(double)*ssd_spec->SEGNUM);
	memset(ssd->VPG, 0., sizeof(double)*ssd_spec->MAXGNUM);
	memset(ssd->EPG, 0.0, sizeof(double)*ssd_spec->MAXGNUM);
	memset(ssd->DEBUG_VPG, 0., sizeof(double)*ssd_spec->MAXGNUM);
	memset(ssd->DEBUG_EPG, 0., sizeof(double)*ssd_spec->MAXGNUM);
	memset(ssd->TMP_VPG, 0.0, sizeof(double)*ssd_spec->MAXGNUM);
	memset(ssd->TMP_EPG, 0.0, sizeof(double)*ssd_spec->MAXGNUM);
	memset(ssd->CUR_VPG, 0.0, sizeof(double)*ssd_spec->MAXGNUM);
	memset(ssd->TOTAL_GNUM, 0.0, sizeof(double)*ssd_spec->MAXGNUM);

	memset(ssd->err_VPG, 0., sizeof(double)*ssd_spec->MAXGNUM);
	memset(ssd->err_EPG, 0.0, sizeof(double)*ssd_spec->MAXGNUM);


	printf("Device size: %.2f GiB\n", (double)(ssd_spec->DEVSIZE/(1024*1024*1024)));
	printf("Logical size: %.2f GB\n", (double)(ssd_spec->LOGSIZE/(1000*1000*1000)));
	printf("# of segments: %d\n", ssd_spec->SEGNUM);
	printf("Segment size: %.2f MB\n", (double)ssd_spec->SEGSIZE/(1024.*1024));
	printf("===========================================================\n");

}

void group_init(struct SSD* ssd, class GROUP *group[], int gnum, int mida_on, int size[]){
	int i;
	int sum = 0;
	//printf("\nGROUP SETUP\n");
	for(i = 0; i < gnum; i++){
		group[i] = new GROUP;
		if (i != gnum-1) {
			//-1: decrease active segment
			group[i]->size = 1; //group size initialization for group 0 ~ n-2
			sum += 1;
			printf("GROUP %d: %d Segments\n", i, group[i]->size);
		}
	}
	group[gnum-1]->size = ssd_spec->SEGNUM - ssd_spec->FREENUM - sum - 1;	// group size init for group n-1
	//consider active segment
	if (mida_on == 1) group[gnum-1]->size -= ssd_spec->NAIVE_N;
	printf("GROUP %d: %d Segments\n", gnum-1, group[gnum-1]->size);

	
	for (int i = 0; i < ssd_spec->GROUPNUM; i++){
		int blk_idx = queue_pop_back(&fbqueue);
		ssd->active[i][0] = blk_idx; //active block per group initializtion
		ssd->active[i][1] = 0; //first active block appending point: 0
		ssd->gnum_info[blk_idx] = i;
	}
}

