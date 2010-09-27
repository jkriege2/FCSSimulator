#ifndef __CORRELATOR_H_INCLUDED__
#define __CORRELATOR_H_INCLUDED__

#include <stdlib.h>
#include <stdio.h>

#include "correlator_block.h"

/**
 * Correlator class implementing a multi-tau correlator.
 * This class is NOT intended to be optimized or fast,
 * it tries to build a hardware equivalent of a parallel
 * multi-tau correlator without beeing parallel. ;-)
 * This is how it should look like (for a (3,2) invocation):
 *
 *
 */
class correlatorjb
{
	private:
		block **blocks;/**Blocks running at different frquencies*/
		unsigned int block_count;/**Number of blocks in current correlator setup*/
		correlatorjb *next;/**One way links list, reference to the next block in chain.*/
	public:
		correlatorjb(unsigned int,unsigned int);
		~correlatorjb();
		void run(int,int);
		void print(unsigned int*);
		unsigned int		get_channel_count();
		int * get_array();
		float * get_array_G();
};

/**
 * Class Constructor, everything is allocted at runtime
 * FIXME: Error handling is still missing!
 * @param block_count Number of blocks to use
 * @param lag_count Number of lags per block
 */
correlatorjb::correlatorjb(unsigned int block_count,unsigned int lag_count)
{
	next=0;
	correlatorjb::block_count=block_count;
	blocks = (block **)malloc(block_count * sizeof(block*));
	if (blocks==0)
		printf("Error!");
	/**Create blocks*/
	for(unsigned int x=0;x<block_count;x++)
		blocks[x]=new block(x,lag_count);
	/**Create linked lists*/
	for(unsigned int x=0;x<(block_count-1);x++)
		blocks[x]->set_next(blocks[x+1]);
}

/**
 * Class Destructor
 */
correlatorjb::~correlatorjb () {
	for(unsigned int x=0;x<block_count;x++)
		delete blocks[x];
}

/**
 * Run the correlator for one step.
 * @param value_global Global value, used for all blocks and lags simultaniousely
 * @param value_local Local value, only used for first block.
 */
void correlatorjb::run(int value_global,int value_local)
{
	if(block_count>0)
		blocks[0]->run(value_global,value_local);
}

/**
 * Print the sums of all lags.
 * TODO: Switch to cout streaming and operator<< overloading
 * @param pos Starting point for printing. Contains the next free position afterwards.
 */
void correlatorjb::print(unsigned int *pos)
{
	for(unsigned int x=0;x<block_count;x++)
		blocks[x]->print(pos);
}

unsigned int correlatorjb::get_channel_count()
{
	unsigned int channel_count=0;
	for(unsigned int x=0;x<block_count;x++)
		channel_count+=blocks[0]->get_lag_count()*(1<<x);
	return channel_count;
}


int * correlatorjb::get_array()
{
	//int channel0=blocks[0]->get_channel(0);
	unsigned int channel_count=this->get_channel_count();
	int *array = (int*)malloc(channel_count * sizeof(int));

	int pos=0;
	for(unsigned int x;x<block_count;x++)
	{
		blocks[x]->fill_array(&(array[pos]));
		pos+=blocks[0]->get_lag_count()*(1<<x);
	}

	return array;
}

float * correlatorjb::get_array_G()
{
	int *array=this->get_array();
	unsigned int channel_count=this->get_channel_count();
	float *result = (float*)malloc(channel_count * sizeof(float));

	float channel0=blocks[0]->get_sum(0);
	for(int x=0;x<(this->get_channel_count());x++)
	{
		result[x]=((float)array[x])/(channel0*channel0);
	}
	delete array;
	return result;
}

#endif
