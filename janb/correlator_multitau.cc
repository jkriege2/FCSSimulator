#include "correlator_multitau.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Class Constructor, everything is allocted at runtime
 * FIXME: Error handling is still missing!
 * @param block_count Number of blocks to use
 * @param lag_count Number of lags per block
 */
correlatorjb::correlatorjb(unsigned int block_count,unsigned int lag_count)
{
	next=0;
	sum=0;
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
    sum+=value_global;
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


int * correlatorjb::get_array_lin()
{
	//int channel0=blocks[0]->get_channel(0);
	unsigned int channel_count=this->get_channel_count();
	int *array = (int*)malloc(channel_count * sizeof(int));

	int pos=0;
	for(unsigned int x;x<block_count;x++)
	{
		blocks[x]->fill_array_lin(&(array[pos]));
		pos+=blocks[0]->get_lag_count()*(1<<x);
	}

	return array;
}

float * correlatorjb::get_array_lin_G()
{
	int *array=this->get_array_lin();
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

int **correlatorjb::get_array()
{
	int **result=(int**)malloc(2 * sizeof(int*));
	result[0]=(int*)malloc(block_count*blocks[0]->get_lag_count() * sizeof(int));;
	result[1]=(int*)malloc(block_count*blocks[0]->get_lag_count() * sizeof(int));;

	unsigned int tau=0;
	for(int x=0;x<block_count;x++)
		blocks[x]->fill_array(result,x*(blocks[0]->get_lag_count()),&tau);

	return result;

}

float **correlatorjb::get_array_G()
{
	int **array=get_array();
	float **result=(float**)malloc(2 * sizeof(float*));
	result[0]=(float*)malloc(block_count*blocks[0]->get_lag_count() * sizeof(float));;
	result[1]=(float*)malloc(block_count*blocks[0]->get_lag_count() * sizeof(float));;
	printf("block_count=%d,    *blocks[0]->get_lag_count()=%d\n", block_count, blocks[0]->get_lag_count());

	float channel0=blocks[0]->get_sum(0);
	for(int x=0;x<(block_count*blocks[0]->get_lag_count());x++)
	{

		result[0][x]=(float)array[0][x];
		result[1][x]=((float)array[1][x])/(sum*sum)*(float)(blocks[x/blocks[0]->get_lag_count()]->get_ticks());
		printf("%d   %ld   %d   %d\n", array[1][x], (long)sum, x/blocks[0]->get_lag_count(), blocks[x/blocks[0]->get_lag_count()]->get_ticks());
	}

	free(array[1]);
	free(array[0]);
	free(array);
	return result;

}

/**
 *Block Constructor. Uses dynamic allocation.
 *FIXME: Error handling is far from perfect.
 *@param id Id of the current block within a correlator. Should range from 0..(x-1).
 *@param lag_count Number of lags in this block.
 */
block::block(unsigned int id,unsigned int lag_count)
{
	next=0;
	ticks=0;
	block::lag_count=lag_count;
	lags = (lag **)malloc(lag_count * sizeof(lag*));
	if (lags==0)
		printf("Error!");
	for(unsigned int x=0;x<lag_count;x++)
		lags[x]=new lag(x,this);
	for(unsigned int x=0;x<(lag_count-1);x++)
		lags[x]->set_next(lags[x+1]);

	block::id=id;
}

/**
 *Class destructor.
 */
block::~block () {
	for(unsigned int x=0;x<lag_count;x++)
		delete lags[x];
}

/**
 * Run the correlation within this block for one time step.
 * @param value_global Global value, used for all blocks and lags simultaniousely
 * @param value_local Local value, only used for first block.
 */
void block::run(int value_global,int value_local)
{
	if(lag_count>0)
		lags[0]->run(value_global,value_local);
	if((ticks%2)==0)
		if(next!=0)
			next->run(value_global,lags[lag_count-1]->acc);
	ticks++;
}

/**
 * Print the sums of all lags.
 * TODO: Switch to cout streaming and operator<< overloading
 * @param pos Starting point for printing. Contains the next free position afterwards.
 */
void block::print(unsigned int *pos)
{
	for(unsigned int x=0;x<lag_count;x++)
	{
		for(int y=0;y<(1 << id);y++)
			lags[x]->print((*pos)++);
	}
}

/**
 *Set the next block in chain. Used for creating linked lists.
 *@param next Next block in chain.
 */
void block::set_next (block *next)
{
	block::next=next;
}

int block::get_lag_count()
{
	return lag_count;
}

void block::fill_array_lin(int *array)
{
	for(int x=0;x<lag_count;x++)
		for(int y=0;y<(1<<id);y++)
			array[x*(1<<id)+y]=lags[x]->get_sum();
}

int block::get_sum(unsigned int lag_number)
{
	if((lag_count>0)&&(lag_number<lag_count))
		return lags[lag_number]->sum;
}

void block::fill_array(int **array,unsigned int pos, unsigned int *tau)
{
	for(unsigned int x=0;x<lag_count;x++)
	{
		array[0][pos+x]=*tau;
		array[1][pos+x]=lags[x]->get_sum();
		*tau+=(1<<id);
	}
}

/**
 * Class Constructor.
 * @param id Id of the current lag in block. Range: 0..(x-1).
 * @param parent Currently unused. Might be used to update acc one klevel above...
 */
lag::lag (unsigned int id,void *parent)
{
	value.old=0;
	value.cur=0;
	next=NULL;
	sum=0;
	acc=0;
	lag::id=id;
}

/**
 *Set the next lag in chain. Used for creating linked lists.
 *@param next Next lag in chain.
 */
void lag::set_next (lag *next)
{
	lag::next=next;
}

/**
 * Run the correlation for this lag for one time step. If this is the last lag in chain, acc is updated.
 * @param value_global Global value, used for all blocks and lags simultaniousely
 * @param value_local Local value, only used for first block.
 */
void lag::run(int value_global, int value_local)
{
	sum+=value_global*value_local;
	if (next!=0)
		next->run(value_global,value.cur);
	else
		acc=(value.cur+value.old);

	/**Implement the FlipFlop*/
	value.old=value.cur;
	value.cur=value_local;
}

/**
 * Print the result of the correlation for this lag.
 * @param pos Lag position in timing context.
 */
void lag::print(unsigned int pos)
{
	printf("%u %i\n",pos,sum);
}

int lag::get_sum()
{
	return sum;
}
