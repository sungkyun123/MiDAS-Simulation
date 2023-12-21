#include <deque>

using namespace std;

int queue_pop_back(deque<int> *blkqueue){
        int blk_idx = blkqueue->back();
        blkqueue->pop_back();
        return blk_idx;
} //pop and identifying popped element in same function

int queue_pop_front(deque<int> *blkqueue){
        int blk_idx = blkqueue->front();
        blkqueue->pop_front();
        return blk_idx;
}


