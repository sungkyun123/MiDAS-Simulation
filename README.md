# MiDAS

## What is MiDAS?
MiDAS is a Migration-based Data placement technique with Adaptive group number and Size configuration for log-structured systems.

MiDAS is implemented in both **trace-driven simulation** and a **real SSD prototype**. The simulator is used to quickly evaluate WAF of the GC techniques. Also, we used a real SSD prototype to measure I/O performance and the overheads associated with the CPU and memory for system execution. You can check the MiDAS's code for the real SSD prototype in the following link (https://github.com/dgist-datalab/FlashFTLDriver/tree/MiDAS). 

The original paper that introduced MiDAS is currently in the revision stage of [USENIX FAST 2024](https://www.usenix.org/conference/fast24).

The archive snapshot with a detailed screencast used in the paper is available at ~~link~~ (not available until Feb 27, 2024).

## Trace-driven simulation
The simulation of MiDAS is implemented by C++.

You can evaluate the WAF of MiDAS with the simulator and trace files.

The MiDAS algorithm is implemented in the directory `model.cpp`.


### Prerequisites
The hardware/software requirements for executing MiDAS are as followed.


#### Hardware
* `DRAM`: Simulation does not require large DRAM capacity. You need 1% of simulated device size for the data structures of MiDAS and the metadata of the device (e.g., L2P mapping table, OOB, invalid/valid bitmap) to test the trace files. For example, you need about 1.3GB size of DRAM to run trace file with 128GB device size.


#### Software
There are special requirements to run the simulation.


### Installation & Compilation
* Clone required reopsitory (MiDAS SSD prototype) into your host machine.
```
$ git clone git@github.com:sungkyun123/MiDAS-Simulation.git
$ cd MiDAS
```

* Download the trace file to test the prototype
```
$ wget https://zenodo.org/record/10409599/files/test-fio-small
```


### Compile & Execution
After downloading trace file, you can test MiDAS. This experiment may be finished in 5 minutes by your server environment. The trace file is suitable for 8GB device.

* Smaller FIO Zipfian 1.0 workload (filename: test-fio-small)
   * You need about 80MB size of DRAM to test this trace file.
```
$ make clean; make
$ ./ssdsimul {trace_file} {victim_selection (fifo/greedy/cost-benefit} {device_size (GB)} {segment_size (MB)}
```

### Some statements for code structure
MiDAS algorithm is implemented in followed files.
This includes UID, MCAM, GCS, and hot block separation algorithms.
- `algorithm.cpp`     : adaptably change group cofiguration, and check irregular pattern
- `model.cpp`     : UID, MCAM, and GCS algorithm
- `hf.cpp`        : Hot block separation


### Results
During the experiment, you can see that MiDAS adaptably change the group configuration.


* UID information : When you run MiDAS, You can see the simulation setup and the parameters of the UID at the beginning.
We sample a subset of LBA for timestamp monitoring, with a sampling rate of 0.01.
We use a coarse-grained update interval unit of 0.25 segments and epoch lengths of 4x of the storage capacity (128GB).
Following result is an example.

```
======================SIMULATION SETUP=====================
===========================================================
Workload: trace/test-fio-small
victim selection: cost-benefit
Device size: 8.00 GiB
Logical size: 8.00 GB
# of segments: 2048
Segment size: 4.00 MB
===========================================================
*** MINIATURE MODEL SETTINGS ***
- lba sampling ratio: 0.01
- interval unit size: 0.250 segments (256 pages)
- time window: 32.00GB (32768 units)
```

* You can see the group configuration and valid ratio of the groups per 1x write request of the storage capacity.

```
[Progress: 72GB] WAF: 2.188, TMP_WAF: 2.392, Utilization: 1.000
HOT[8]: 0.0201 (ERASE: 1359)  (Q size: 8) (desig_size: 8)
GROUP 1[49]: 0.5334 (ERASE: 716)  (Q size: 49) (desig_size: 49)
GROUP 2[79]: 0.7010 (ERASE: 382)  (Q size: 79) (desig_size: 79)
GROUP 3[93]: 0.7691 (ERASE: 268)  (Q size: 93) (desig_size: 93)
GROUP 4[97]: 0.8040 (ERASE: 185)  (Q size: 97) (desig_size: 97)
--------------------------
GROUP 5[1712]: 0.9146 (ERASE: 1988)  (Q size: 1712) (desig_size: 1712)

```


* When an epoch is over, GCS algorithm finds the best group configuration using UID and MCAM. The group configuration is shown as follows.

```
*****MODEL PREDICTION RESULTS*****
group number: 6
*group 0: size 9, valid ratio 0.000000
*group 1: size 51, valid ratio 0.560998
*group 2: size 80, valid ratio 0.706775
*group 3: size 94, valid ratio 0.780401
*group 4: size 98, valid ratio 0.838156
*group 5: size 1712, valid ratio 0.908237
calculated WAF: 2.220729
Used traffic (CALC) : 0.702
************************************
```


* MiDAS periodically check the irregular pattern of the workload. If there is a group that its valid ratio prediction is wrong, MiDAS gives up on adjusting group sizes for all groups beyond the group and simply merges the groups.

```
**Error comparing function**
[group 1]'s real vs exp valid ratio:[0.5333 vs 0.5610]
[group 2]'s real vs exp valid ratio:[0.7008 vs 0.7068]
[group 3]'s real vs exp valid ratio:[0.7697 vs 0.7804]
[group 4]'s real vs exp valid ratio:[0.8060 vs 0.8382]
[group 5]'s real vs exp valid ratio:[0.9141 vs 0.9082]
```


* When the excution is over, You can check the total user-write, GC-write, WAF, and the average valid data ratio of each group on the MiDAS simulation.
```
 Experimental result
============================================
USERWRITE: 54237899
COPYWRITE: 70485894
TOTAL WAF: 2.300
HOT[9]: 0.0202 (ERASE: 22858)
GROUP 1[47]: 0.5408 (ERASE: 12136)
GROUP 2[73]: 0.7062 (ERASE: 6559)
GROUP 3[87]: 0.7721 (ERASE: 4632)
GROUP 4[91]: 0.8046 (ERASE: 3572)
GROUP 5[1737]: 0.9090 (ERASE: 31715)
```
